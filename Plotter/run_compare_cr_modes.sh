#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
# ENV=("StandAlone" "CMSSW")
ENV=("StandAlone")
SYSTEMS=("pp")
ENERGIES=("13TeV")
PTHATS=("10")
# LDMES=("1")
LDMES=("100")
# EVENTS=("10000000")
EVENTS=("200000")
RIVET_CONFIGS=("CMS_default" "LHCb_default" "CMS_FineBin_R3" "CMS_FineBin_R4" "CMS_20Bin_R3" "CMS_20Bin_R4" "CMS_50Bin_R3" "CMS_50Bin_R4" "CMS_100Bin_R3" "CMS_100Bin_R4")

for ENV in "${ENV[@]}"; do
  for SYSTEM in "${SYSTEMS[@]}"; do
    for ENERGY in "${ENERGIES[@]}"; do
      for PTHAT in "${PTHATS[@]}"; do
        for LDME in "${LDMES[@]}"; do
          for EVENT in "${EVENTS[@]}"; do
            for RIVET_CONFIG in "${RIVET_CONFIGS[@]}"; do
              OUTPUT_FOLDER="${SCRIPT_DIR}/fig/CRMode_Comparisons/${ENV}_${SYSTEM}_${ENERGY}_LDMEfac${LDME}_ptHat${PTHAT}_${EVENT}_rivet_config_${RIVET_CONFIG}"
              python3 "${SCRIPT_DIR}/compare_cr_modes.py" \
                --energy "${ENERGY}" \
                --system "${SYSTEM}" \
                --ptHat "${PTHAT}" \
                --ldme "${LDME}" \
                --event "${EVENT}" \
                --rivet-config "${RIVET_CONFIG}" \
                --input-template "${ROOT_DIR}/Production/rivetFile/${ENV}/{Energy_value}/{System_value}_OniaShower{channel}_{Energy_value}_ptHat{ptHat_value}_LDMEfac{LDME_value}_CR{CRMode}_{Event_value}.yoda" \
                --output-folder "${OUTPUT_FOLDER}" \
                "$@"
            done
          done
        done
      done
    done
  done
done
