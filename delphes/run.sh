#!/bin/bash

set -euo pipefail

if [ -z "${4:-}" ]; then
    echo "Must enter 4 arguments"
    echo -e "\t1: Dataset Tag (from Pythia/HepMC input)"
    echo -e "\t2: Dataset Tag (for Delphes ROOT output)"
    echo -e "\t3: Number of runs"
    echo -e "\t4: Max num cpu cores"
    exit 1
fi

PY_tag=$1
DELPHES_tag=$2
runs=$3
max_cpu_cores=$4

small_r_algorithm=${DELPHES_SMALL_R_ALGORITHM:-antikt}
small_r=${DELPHES_SMALL_R:-0.4}
small_r_ptmin=${DELPHES_SMALL_R_PTMIN:-10.0}
large_r_algorithm=${DELPHES_LARGE_R_ALGORITHM:-cambridge}
large_r=${DELPHES_LARGE_R:-1.5}
large_r_ptmin=${DELPHES_LARGE_R_PTMIN:-250.0}
small_r_jet_flavor_dr=${DELPHES_SMALL_R_JET_FLAVOR_DR:-0.4}
small_r_parton_ptmin=${DELPHES_SMALL_R_PARTON_PTMIN:-1.0}
small_r_parton_eta_max=${DELPHES_SMALL_R_PARTON_ETA_MAX:-2.5}
card_template=${DELPHES_CARD_TEMPLATE:-delphes_card_ATLAS_top_polarimetry.tcl.template}

workspace="WS_${DELPHES_tag}"
mkdir -p "$workspace"
mkdir -p "$workspace/logs"
mkdir -p "$workspace/config"

rendered_card="$workspace/config/delphes_card_${DELPHES_tag}.tcl"
python render_card.py \
    --template "$card_template" \
    --output "$rendered_card" \
    --small-r-algorithm "$small_r_algorithm" \
    --small-r "$small_r" \
    --small-r-ptmin "$small_r_ptmin" \
    --large-r-algorithm "$large_r_algorithm" \
    --large-r "$large_r" \
    --large-r-ptmin "$large_r_ptmin" \
    --small-r-jet-flavor-dr "$small_r_jet_flavor_dr" \
    --small-r-parton-ptmin "$small_r_parton_ptmin" \
    --small-r-parton-eta-max "$small_r_parton_eta_max"

echo "Please be patient while Delphes processes the HepMC files..."

job=0
batch=1
pids=()
for (( i=0 ; i<$runs ; i++ ));
do
    input_file="../pythia/WS_${PY_tag}/dataset_${PY_tag}_${i}.hepmc"
    output_file="$workspace/dataset_${DELPHES_tag}_${i}.root"
    log_file="$workspace/logs/dataset_${DELPHES_tag}_${i}.log"

    echo -e "\t\tSubmitting job to Delphes run $i"
    (
        set -euo pipefail
        trap 'rm -f "$output_file"' ERR

        if [[ ! -s "$input_file" ]]; then
            echo "Missing or empty HepMC input: $input_file" >&2
            exit 1
        fi

        ../submodules/delphes-v3.5.0/DelphesHepMC2 \
            "$rendered_card" \
            "$output_file" \
            "$input_file"

        python validate_delphes_output.py --hepmc "$input_file" --root "$output_file"
    ) > "$log_file" 2>&1 &
    pids+=("$!")

    job=$((job+1))
    if [ "$job" == "$max_cpu_cores" ]; then
        echo -e "\tStopping job submissions! Please wait for batch $batch to finish..."
        for pid in "${pids[@]}"; do
            wait "$pid" || exit 1
        done
        pids=()
        job=0
        batch=$((batch+1))
    fi
done

for pid in "${pids[@]}"; do
    wait "$pid" || exit 1
done

echo -e "\tDelphes Simulation Done!"
