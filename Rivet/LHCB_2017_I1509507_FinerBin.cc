// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"


namespace Rivet {


  /// @brief J/psi production in jets at 13 TeV
  class LHCB_2017_I1509507_FinerBin : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2017_I1509507_FinerBin);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      declare(FinalState(), "FS");
      // Define the range and binning parameters for the fragmentation variable z
      const double Z_MIN = 0.0;
      const double Z_MAX = 1.0;
      const int NUM_BINS = 100; // For a bin width of (1.0 - 0.0) / 100 = 0.01
      
      // Book histograms with custom binning
      book(_h_frag[0], "frag_prompt",    NUM_BINS, Z_MIN, Z_MAX);
      book(_h_frag[1], "frag_nonprompt", NUM_BINS, Z_MIN, Z_MAX);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // first see if we have any J/psi in the region
      Particles Jpsi;
      for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==443 &&
										Cuts::etaIn(2.,4.5))) {
	if(p.children().size()!=2) continue;
	bool found=true;
	for(const Particle & child :p.children()) {
	  double eta = child.eta();
	  if (child.abspid()!=13 || eta<2. || eta>4.5 || child.perp()<0.5 || child.p3().mod()<5.) {
	    found=false;
	    break;
	  }
	}
	if (found) Jpsi.push_back(p);
      }
      // no jpsi veto
      if(Jpsi.empty()) vetoEvent;
      // now get the final-state particles for the jets
      const Particles & fsParticles = apply<FinalState>(event, "FS").particles();
      vector<PseudoJet> particles;
      particles.reserve(fsParticles.size());
      // fs for fastjet omitting any J/psi decay products
      for(const Particle & p : fsParticles) {
	// skip muons and neutrinos (seems standard for LHCb)
	if(p.abspid()==13 || p.abspid()==12 || p.abspid()==14 or p.abspid()==16) continue;
	// skip anything coming from the decay of one of the jpsis
	Particle parent = p;
	while(!parent.parents().empty()) {
	  if(parent.pid()==443) break;
	  parent=parent.parents()[0];
	}
	bool match = parent.pid()==443;
	if(match) {
	  match =false;
	  for(const Particle & psi : Jpsi) {
	    match = fuzzyEquals(parent.momentum(),psi.momentum());
	    if(match) break;
	  }
	}
	if(!match) {
	  PseudoJet j = p.pseudojet();
	  j.set_user_index(0);
	  particles.push_back(j);
	}
      }
      // add the jpsis to the particles for fastjet
      for(const Particle & p : Jpsi) {
	PseudoJet j = p.pseudojet();
	j.set_user_index(p.fromBottom()+1);
	particles.push_back(j);
      }
      JetDefinition jet_def(fastjet::antikt_algorithm, 0.5);
      fastjet::ClusterSequence clu = ClusterSequence(particles,jet_def); 
      vector<PseudoJet> jets = clu.inclusive_jets();
      for(const PseudoJet & jet : jets) {
	// pt and eta cut
	if(jet.perp()<20. || jet.eta()<2. || jet.eta()>4.) continue;
	// loop over constituents and find jpsi
	for(const PseudoJet & sub : jet.constituents()) {
	  if(sub.user_index()==0) continue;
	  double z = sub.perp()/jet.perp();
	  _h_frag[sub.user_index()-1]->fill(z);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h_frag[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_frag[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2017_I1509507_FinerBin);

}
