#!/bin/bash

# Read options from config
if [[ -z $1 ]]; then
  source job.config
else
  source $1
fi

# Get current directory
WORKING_DIR=$(pwd)

# Setup python venv
source setup.sh

# Run MadGraph Generation
if [ "$bypass_madgraph" = false ]; then
  start=`date +%s`
  cd madgraph
  ./run.sh "$MG_tag" "$process" "$generation" "$num_runs" "$num_events_per_run" "$max_cpu_cores" "$seed" "$top_pT_cut"
  cd $WORKING_DIR
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\tTime (sec): $runtime"
fi

# Calculate number of batches based on max number of cpu available
# Use python to perform calculations
num_batches=$(python -c "print(int($num_runs/float($max_cpu_cores))) if $num_runs%$max_cpu_cores==0 else print(int($num_runs/float($max_cpu_cores))+1)")

# Run Pythia Showering
if [ "$bypass_pythia" = false ]; then
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
    ./run.sh "$PY_tag" "$DELPHES_tag" "$num_runs" "$max_cpu_cores"
  )
  delphes_status=$?
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\tDelphes Time (sec): $runtime"
  if [ "$delphes_status" -ne 0 ]; then
    exit "$delphes_status"
  fi
fi

# Run preprocessing to match fat jets to top quarks
if [ "$bypass_preprocessing" = false ]; then
  echo "Please be patient for Preprocessing..."
  start=`date +%s`
  cd model

  job=0
  batch=1
  for (( i=0 ; i<$num_runs ; i++ ));
  do
    echo -e "\t\tSubmitting job to preprocess run $i"
    mkdir -p "${dir_datasets}/run_$i"
    python -u DataLoader_Parallel.py $PY_tag $i $dir_datasets > "${dir_datasets}/run_$i/preprocessing.log" &
    job=$((job+1))
    if [ $job == $max_cpu_cores ]; then
      echo -e "\tStopping jobs submissions! Please wait for batch $batch to finish..."
      wait
      job=0
      batch=$((batch+1))
    fi
  done
  wait

  python -u Combine_DataSets.py $PY_tag $num_runs $dir_datasets
  cd $WORKING_DIR
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\tBatching Done!"
  echo -e "\tTime (sec): $runtime"
fi

# Run training script
if [ "$bypass_train" = false ]; then
  echo "Please be patient for Training..."
  start=`date +%s`
  cd model
  mkdir -p $dir_training
  mkdir -p "$dir_training/models"
  python -u New_Training.py $PY_tag $epochs $embed_dim $dir_datasets $dir_training $analysis_type | tee "${dir_training}/training.log"
  cd $WORKING_DIR
  end=`date +%s`
  runtime=$((end-start))
  echo -e "\tTraining Done!"
  echo -e "\tTime (sec): $runtime"
fi
