// main132_rivet.cc
// Direct Pythia + Rivet driver that writes YODA without an intermediate HepMC file.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.

#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/InputParser.h"
#include "Pythia8Plugins/HepMC3.h"

#include "Rivet/AnalysisHandler.hh"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <mutex>
#include <atomic>
#include <string>
#include <vector>

using namespace Pythia8;

namespace {

// Protect Rivet analysis and YODA writing from concurrent access.
std::mutex rivetMutex;
// Timing counters to diagnose lock contention (nanoseconds)
static std::atomic<long long> rivet_lock_wait_ns{0};
static std::atomic<long long> rivet_locked_ns{0};

long countJpsi(const Pythia8::Event& event) {
  long count = 0;
  for (const auto& particle : event) {
    if (particle.idAbs() == 443) ++count;
  }
  return count;
}

const std::vector<std::string> analyses = {
  "LHCB_2017_I1509507",
  "LHCB_2017_I1509507_FinerBin",
  "LHCB_2017_I1509507_CMSCut_R4",
  "LHCB_2017_I1509507_CMSCut_R3",
  "CMS_2022_I1870319",
  "CMS_2022_I1870319_FineBin_R3",
  "CMS_2022_I1870319_FineBin_R4",
  "CMS_2022_I1870319_20Bin_R3",
  "CMS_2022_I1870319_20Bin_R4",
  "CMS_2022_I1870319_50Bin_R3",
  "CMS_2022_I1870319_50Bin_R4",
  "CMS_2022_I1870319_100Bin_R3",
  "CMS_2022_I1870319_100Bin_R4",
  "CMS_2022_I1870319_LHCbcut_FineBin"
};

} // namespace

//=========================================================================
int main(int argc, char* argv[]) {

  int nevts = -1;
  int nthreads = -1;
  double ldmeFac = -1.0;
  double pthatmin = -1.0;

  std::vector<char*> filtered_argv;
  filtered_argv.push_back(argv[0]);
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "--nevts" || arg == "--nthreads" || arg == "--ldmeFac" || arg == "--pthatmin") && i + 1 < argc) {
      if (arg == "--nevts") nevts = std::stoi(argv[++i]);
      else if (arg == "--nthreads") nthreads = std::stoi(argv[++i]);
      else if (arg == "--ldmeFac") ldmeFac = std::stod(argv[++i]);
      else if (arg == "--pthatmin") pthatmin = std::stod(argv[++i]);
    } else {
      filtered_argv.push_back(argv[i]);
    }
  }

  int filtered_argc = static_cast<int>(filtered_argv.size());

  InputParser ip(
    "This program illustrates how Pythia8 events can be converted directly into Rivet YODA output.",
    {"./main132_rivet -c main132.cmnd -o main132.yoda [--nevts N] [--nthreads N] [--ldmeFac X] [--pthatmin X]"}
  );
  ip.require("c", "Use this user-written command file.", {"-cmnd"});
  ip.require("o", "Specify YODA output filename.", {"-out"});

  const InputParser::Status status = ip.init(filtered_argc, filtered_argv.data());
  if (status != InputParser::Valid) return status;

  const std::string cmnd = ip.get<std::string>("c");
  const std::string out = ip.get<std::string>("o");

  std::cout << "\n >>> PYTHIA settings will be read from file '" << cmnd
            << "' <<< \n >>> Rivet YODA output will be written to file '"
            << out << "' <<< \n";

  PythiaParallel pythia;
  pythia.readFile(cmnd);

  std::cout << "[CMND] Parallelism:numThreads = " << pythia.settings.mode("Parallelism:numThreads") << std::endl;
  std::cout << "[CMND] PhaseSpace:pTHatMin = " << pythia.settings.parm("PhaseSpace:pTHatMin") << std::endl;
  std::cout << "[CMND] Main:numberOfEvents = " << pythia.settings.mode("Main:numberOfEvents") << std::endl;
  std::cout << "[CMND] OniaShower:ldmeFac = " << pythia.settings.parm("OniaShower:ldmeFac") << std::endl;

  if (nthreads > 0) {
    pythia.readString("Parallelism:numThreads = " + std::to_string(nthreads));
    std::cout << "[Override] Parallelism:numThreads = " << nthreads << std::endl;
  }
  if (pthatmin > 0.0) {
    pythia.readString("PhaseSpace:pTHatMin = " + std::to_string(pthatmin));
    std::cout << "[Override] PhaseSpace:pTHatMin = " << pthatmin << std::endl;
  }
  if (nevts > 0) {
    pythia.readString("Main:numberOfEvents = " + std::to_string(nevts));
    std::cout << "[Override] Main:numberOfEvents = " << nevts << std::endl;
  }
  if (ldmeFac > 0.0) {
    pythia.readString("OniaShower:ldmeFac = " + std::to_string(ldmeFac));
    std::cout << "[Override] OniaShower:ldmeFac = " << ldmeFac << std::endl;
  }

  const int nEvent = pythia.settings.mode("Main:numberOfEvents");
  if (!pythia.init()) return 1;

  Pythia8ToHepMC converter;
  Rivet::AnalysisHandler rivet;
  rivet.setCheckBeams(false);
  rivet.addAnalyses(analyses);

  bool rivetInitialized = false;
  long nJpsiProduced = 0;

  const auto runStart = std::chrono::steady_clock::now();

  std::vector<long> eventsPerThread = pythia.run(nEvent, [&](Pythia* pythiaPtr) {
    // measure time waiting for the rivet mutex
    const auto t_before_lock = std::chrono::steady_clock::now();
    std::unique_lock<std::mutex> lock(rivetMutex);
    const auto t_locked = std::chrono::steady_clock::now();
    rivet_lock_wait_ns.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t_locked - t_before_lock).count(), std::memory_order_relaxed);

    // perform the Rivet-protected work and time it
    const auto t_work_start = std::chrono::steady_clock::now();
    nJpsiProduced += countJpsi(pythiaPtr->event);
    if (!converter.fillNextEvent(*pythiaPtr)) {
      std::cerr << "[ERROR] Failed to convert the current Pythia event for Rivet analysis." << std::endl;
      lock.unlock();
      return;
    }
    if (!rivetInitialized) {
      rivet.init(converter.event());
      rivetInitialized = true;
    }
    rivet.analyze(converter.event());
    const auto t_work_end = std::chrono::steady_clock::now();
    rivet_locked_ns.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(t_work_end - t_work_start).count(), std::memory_order_relaxed);
    // lock is released on scope exit
  });

  if (!rivetInitialized) {
    std::cerr << "[ERROR] No events were analyzed, so no YODA file was written." << std::endl;
    return 1;
  }

  const auto runEnd = std::chrono::steady_clock::now();
  const double elapsedSeconds = std::chrono::duration<double>(runEnd - runStart).count();
  const long nGenerated = std::accumulate(eventsPerThread.begin(), eventsPerThread.end(), 0L);
  long nTried = 0;
  long nAccepted = 0;
  double sigmaErr = 0.0;
  bool sigmaErrSet = false;

  pythia.foreach([&](Pythia* pythiaPtr) {
    nTried += pythiaPtr->info.nTried();
    nAccepted += pythiaPtr->info.nAccepted();
    if (!sigmaErrSet) {
      sigmaErr = pythiaPtr->info.sigmaErr();
      sigmaErrSet = true;
    }
  });

  const double efficiency = (nTried > 0) ? static_cast<double>(nAccepted) / static_cast<double>(nTried) : 0.0;
  const double productionRate = (elapsedSeconds > 0.0) ? static_cast<double>(nAccepted) / elapsedSeconds : 0.0;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "[SUMMARY] Events generated = " << nGenerated << std::endl;
  std::cout << "[SUMMARY] J/psi produced = " << nJpsiProduced << std::endl;
  std::cout << "[SUMMARY] Events tried = " << nTried << std::endl;
  std::cout << "[SUMMARY] Events accepted = " << nAccepted << std::endl;
  std::cout << "[SUMMARY] Efficiency = " << efficiency << std::endl;
  std::cout << "[SUMMARY] Production rate = " << productionRate << " events/s" << std::endl;
  std::cout << "[SUMMARY] Cross section = " << pythia.sigmaGen();
  if (sigmaErrSet) {
    std::cout << " +/- " << sigmaErr;
  }
  std::cout << " mb" << std::endl;
  std::cout << "[SUMMARY] Wall time = " << elapsedSeconds << " s" << std::endl;

  // Print diagnostic timing for rivet mutex contention
  {
    const long long wait_ns = rivet_lock_wait_ns.load(std::memory_order_relaxed);
    const long long locked_ns = rivet_locked_ns.load(std::memory_order_relaxed);
    const double wait_s = static_cast<double>(wait_ns) / 1e9;
    const double locked_s = static_cast<double>(locked_ns) / 1e9;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[DEBUG] Rivet lock total wait = " << wait_s << " s, total locked = " << locked_s << " s" << std::endl;
    if (nGenerated > 0) {
      std::cout << "[DEBUG] per-generated-event: wait = " << (wait_s / nGenerated) << " s, locked = " << (locked_s / nGenerated) << " s" << std::endl;
    }
  }

  rivet.finalize();
  rivet.writeData(out);
  pythia.stat();
  return 0;
}