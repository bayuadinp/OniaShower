#!/bin/bash

set -euo pipefail

# --- CONFIGURATION SECTION ---
# Direct YODA production with the in-memory Pythia + Rivet driver.

System=("pp")
# CHANNELS=("Swave" "Pwave" "Singlet" "Octet")
# PTHATS=("20" "25" "10")
CHANNELS=("Singlet")
PTHATS=("25")
CR=("CR0")
Energy=("13TeV")
# Num=("10000000")
Num=("2000000")
LDME=("100")

# 2. Set environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate rivet-mpi

# Remove conda lib from LD_LIBRARY_PATH for system compatibility
export LD_LIBRARY_PATH=$(printf '%s' "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | awk -v remove="${CONDA_PREFIX}/lib" 'length($0) && $0 != remove { print }' | paste -sd: -)

# Add Pythia8 library path for main132_rivet
export LD_LIBRARY_PATH=/home/bayuadinp/software/pythia8316/lib:$LD_LIBRARY_PATH

python --version

set +u
source /home/bayuadinp/software/rivetbootstrap/Rivet-4.1.2/rivetenv.sh
set -u
export RIVET_ANALYSIS_PATH=/home/bayuadinp/software/rivetbootstrap/Rivet-4.1.2/analyses/pluginLHCb/:/home/bayuadinp/software/rivetbootstrap/Rivet-4.1.2/analyses/pluginCMS/:${RIVET_ANALYSIS_PATH:-}

# 3. Ensure output directories exist
mkdir -p rivetFile logs reports

# --- EXECUTION LOOP ---

for CHANNEL in "${CHANNELS[@]}"; do
 for e in "${Energy[@]}"; do
  for PTHAT in "${PTHATS[@]}"; do
   for cr in "${CR[@]}"; do
    for ldme in "${LDME[@]}"; do
     for num in "${Num[@]}"; do
        BASENAME="${System}_OniaShower${CHANNEL}_${e}_ptHat${PTHAT}_LDMEfac${ldme}_${cr}_${num}"
        CMND="${System}_OniaShower${CHANNEL}_${e}_ptHat${PTHAT}_${cr}"
        LOG_FILE="logs/${BASENAME}.log"
        REPORT_FILE="reports/${BASENAME}.report"
        rm -f "$LOG_FILE"
        rm -f "$REPORT_FILE"

        echo "========================================================"
        echo "Processing: $BASENAME $(date)" | tee -a "$LOG_FILE"
        echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
        echo "========================================================"

        if [ -f "cmndFile/${CMND}.cmnd" ]; then
            mkdir -p "rivetFile/StandAlone/${e}"
            echo "--> Running direct Pythia + Rivet generation..." | tee -a "$LOG_FILE"
            nice -n 19 ./main132_rivet \
               -c "cmndFile/${CMND}.cmnd" \
               -o "rivetFile/StandAlone/${e}/${BASENAME}.yoda" \
               --nevts "${num}" --nthreads 1 --ldmeFac "${ldme}" --pthatmin "${PTHAT}" 2>&1 | tee -a "$LOG_FILE"
            if [ ${PIPESTATUS[0]} -ne 0 ]; then
                echo "[ERROR] Direct Rivet production failed for $BASENAME. Check log: $LOG_FILE" | tee -a "$LOG_FILE"
            fi
        else
            echo "[ERROR] Input file cmndFile/${CMND}.cmnd not found!" | tee -a "$LOG_FILE"
            continue
        fi

        {
            echo "========================================================"
            echo "Summary report for: $BASENAME"
            echo "Source log: $LOG_FILE"
            if ! grep '^\[SUMMARY\]' "$LOG_FILE"; then
                echo "[WARN] No summary block found in $LOG_FILE"
            fi
            echo "========================================================"
            echo ""
        } >> "$REPORT_FILE"

        echo "Done with $BASENAME" | tee -a "$LOG_FILE"
        echo "Output YODA file: rivetFile/StandAlone/${e}/${BASENAME}.yoda" | tee -a "$LOG_FILE"
        echo "log file: $LOG_FILE" | tee -a "$LOG_FILE"
        echo "report file: $REPORT_FILE" | tee -a "$LOG_FILE"
        echo ""
     done
    done
   done
  done
 done
done

echo "All jobs completed successfully."