#!/bin/bash

# Read options from config
if [[ -z $1 ]]; then
  source job.config || exit $?
else
  source "$1" || exit $?
fi

stream_hepmc=${stream_hepmc:-false}
case "$stream_hepmc" in
  true|false) ;;
  *) echo "stream_hepmc must be true or false" >&2; exit 1 ;;
esac
if [[ "$stream_hepmc" == true && "$bypass_pythia" != "$bypass_delphes" ]]; then
  echo "With stream_hepmc=true, Pythia and Delphes must both run or both be bypassed." >&2
  exit 1
fi

# Setup python venv
source setup.sh || exit $?

# Run MadGraph Generation
if [ "$bypass_madgraph" = false ]; then
  start=`date +%s`
  (
    cd madgraph || exit 1
    ./run.sh "$MG_tag" "$process" "$generation" "$num_runs" "$num_events_per_run" "$max_cpu_cores" "$seed" "$top_pT_cut"
  ) || exit $?
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\tTime (sec): $runtime"
fi

# Run Pythia separately only in file-based mode.
if [[ "$bypass_pythia" == false && "$stream_hepmc" == false ]]; then
  start=`date +%s`
  (
    cd pythia || exit 1
    ./run.sh "$MG_tag" "$PY_tag" "$num_runs" "$max_cpu_cores"
  )
  pythia_status=$?
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\tPythia Time (sec): $runtime"
  if [ "$pythia_status" -ne 0 ]; then
    exit "$pythia_status"
  fi
fi

# Run Delphes detector simulation
if [ "$bypass_delphes" = false ]; then
  : "${DELPHES_tag:?Missing DELPHES_tag in config}"
  : "${delphes_card_template:?Missing delphes_card_template in config}"
  : "${delphes_small_r_algorithm:?Missing delphes_small_r_algorithm in config}"
  : "${delphes_small_r:?Missing delphes_small_r in config}"
  : "${delphes_small_r_jet_ptmin:?Missing delphes_small_r_jet_ptmin in config}"
  : "${delphes_large_r_algorithm:?Missing delphes_large_r_algorithm in config}"
  : "${delphes_large_r:?Missing delphes_large_r in config}"
  : "${delphes_large_r_jet_ptmin:?Missing delphes_large_r_jet_ptmin in config}"
  : "${delphes_small_r_jet_flavor_dr:?Missing delphes_small_r_jet_flavor_dr in config}"
  : "${delphes_small_r_parton_ptmin:?Missing delphes_small_r_parton_ptmin in config}"
  : "${delphes_small_r_parton_eta_max:?Missing delphes_small_r_parton_eta_max in config}"
  start=`date +%s`
  delphes_args=("$PY_tag" "$DELPHES_tag" "$num_runs" "$max_cpu_cores")
  stage="Delphes"
  if [[ "$stream_hepmc" == true ]]; then
    delphes_args+=("${MG_tag:?Missing MG_tag in config}")
    stage="Pythia + Delphes"
  fi
  (
    cd delphes || exit 1
    DELPHES_CARD_TEMPLATE="$delphes_card_template" \
    DELPHES_SMALL_R_ALGORITHM="$delphes_small_r_algorithm" \
    DELPHES_SMALL_R="$delphes_small_r" \
    DELPHES_SMALL_R_PTMIN="$delphes_small_r_jet_ptmin" \
    DELPHES_LARGE_R_ALGORITHM="$delphes_large_r_algorithm" \
    DELPHES_LARGE_R="$delphes_large_r" \
    DELPHES_LARGE_R_PTMIN="$delphes_large_r_jet_ptmin" \
    DELPHES_SMALL_R_JET_FLAVOR_DR="$delphes_small_r_jet_flavor_dr" \
    DELPHES_SMALL_R_PARTON_PTMIN="$delphes_small_r_parton_ptmin" \
    DELPHES_SMALL_R_PARTON_ETA_MAX="$delphes_small_r_parton_eta_max" \
    ./run.sh "${delphes_args[@]}"
  )
  delphes_status=$?
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\t$stage Time (sec): $runtime"
  if [ "$delphes_status" -ne 0 ]; then
    exit "$delphes_status"
  fi
fi

# Join corresponding Delphes and MadGraph entries and select reconstructed events.
if [ "${bypass_analysis:-true}" = false ]; then
  : "${DELPHES_tag:?Missing DELPHES_tag in config}"
  start=`date +%s`
  ELECTRON_PT_MIN="${ELECTRON_PT_MIN:-30}" \
  ELECTRON_ABS_ETA_MAX="${ELECTRON_ABS_ETA_MAX:-3}" \
  MUON_PT_MIN="${MUON_PT_MIN:-30}" \
  MUON_ABS_ETA_MAX="${MUON_ABS_ETA_MAX:-3}" \
  MET_MIN="${MET_MIN:-30}" \
  LARGE_R_JET_PT_MIN="${LARGE_R_JET_PT_MIN:-250}" \
  LARGE_R_JET_ABS_ETA_MAX="${LARGE_R_JET_ABS_ETA_MAX:-3}" \
  JET_OVERLAP_DR="${JET_OVERLAP_DR:-1.5}" \
  LEPTONIC_B_JET_PT_MIN="${LEPTONIC_B_JET_PT_MIN:-30}" \
  LEPTONIC_B_JET_ABS_ETA_MAX="${LEPTONIC_B_JET_ABS_ETA_MAX:-3}" \
  LEPTON_B_JET_DR_MIN="${LEPTON_B_JET_DR_MIN:-0.3}" \
  LEPTON_B_JET_DR_MAX="${LEPTON_B_JET_DR_MAX:-3.0}" \
  MAX_LARGE_R_CONSTITUENTS="${MAX_LARGE_R_CONSTITUENTS:-100}" \
  MAX_SMALL_R_JETS="${MAX_SMALL_R_JETS:-10}" \
  bash analysis/run.sh "$MG_tag" "$DELPHES_tag" "$num_runs" "$max_cpu_cores" || exit $?
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\tAnalysis Time (sec): $runtime"
fi
