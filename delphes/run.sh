#!/bin/bash

set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 || ! $3 =~ ^[1-9][0-9]*$ || ! $4 =~ ^[1-9][0-9]*$ ]]; then
    echo "Must enter 4 arguments (and optionally a MadGraph tag to stream)"
    echo -e "\t1: Dataset Tag (from Pythia/HepMC input)"
    echo -e "\t2: Dataset Tag (for Delphes ROOT output)"
    echo -e "\t3: Number of runs"
    echo -e "\t4: Max num cpu cores"
    echo -e "\t5: MadGraph tag (optional; stream Pythia into Delphes)"
    exit 1
fi

PY_tag=$1
DELPHES_tag=$2
runs=$3
max_cpu_cores=$4
stream_mg_tag=${5:-}

if [[ -n "$stream_mg_tag" ]] && (( runs > max_cpu_cores )); then
    echo "Streaming requires num_runs <= max_cpu_cores (got $runs > $max_cpu_cores)." >&2
    exit 1
fi

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

if [[ -n "$stream_mg_tag" ]]; then
    make -C ../pythia/src generate_dataset
    echo "Streaming Pythia into Delphes ($runs concurrent runs)..."
else
    echo "Please be patient while Delphes processes the HepMC files..."
fi

delphes_bin=../submodules/delphes-v3.5.0/DelphesHepMC2

stream_run() (
    local run=$1 output_file=$2 log_file=$3
    local label_file="../madgraph/pp_tt_semi_full_${stream_mg_tag}/labels_${stream_mg_tag}_${run}.root"
    local pythia_log="$workspace/logs/pythia_${DELPHES_tag}_${run}.log"
    local tmp_dir tmp_root statuses

    if [[ ! -s "$label_file" ]]; then
        echo "Missing or empty MadGraph labels for run $run: $label_file" >&2
        exit 1
    fi

    tmp_dir=$(mktemp -d "$workspace/.stream_${run}.XXXXXX")
    tmp_root="$tmp_dir/dataset.root"  # Delphes requires a nonexistent output file.
    trap 'rm -f -- "$tmp_root"; rmdir -- "$tmp_dir"' EXIT

    # Keep HepMC alone in the pipe: fd 3 saves Pythia's piped stdout before
    # stdout and stderr go to its log. Delphes's '-' reads stdin without seeking.
    if (
        cd ../pythia/src || exit 1
        exec ./run_dataset "$stream_mg_tag" "$PY_tag" "$run" /proc/self/fd/3
    ) 3>&1 >"$pythia_log" 2>&1 |
        "$delphes_bin" "$rendered_card" "$tmp_root" - >"$log_file" 2>&1
    then
        python validate_delphes_output.py --pythia-log "$pythia_log" --labels "$label_file" \
            --root "$tmp_root" >>"$log_file" 2>&1
        mv -fT -- "$tmp_root" "$output_file"
    else
        statuses=("${PIPESTATUS[@]}")
        echo "Streaming run $run failed (Pythia=${statuses[0]}, Delphes=${statuses[1]}). See $pythia_log and $log_file." >&2
        exit 1
    fi
)

wait_workers() {
    local i failed=0
    for i in "${!pids[@]}"; do
        if ! wait "${pids[i]}"; then
            echo "Delphes job failed; see ${logs[i]}" >&2
            failed=1
        fi
    done
    pids=()
    logs=()
    return "$failed"
}

job=0
batch=1
pids=()
logs=()
for (( i=0 ; i<$runs ; i++ ));
do
    input_file="../pythia/WS_${PY_tag}/dataset_${PY_tag}_${i}.hepmc"
    output_file="$workspace/dataset_${DELPHES_tag}_${i}.root"
    log_file="$workspace/logs/dataset_${DELPHES_tag}_${i}.log"

    echo -e "\t\tSubmitting job to Delphes run $i"
    if [[ -n "$stream_mg_tag" ]]; then
        stream_run "$i" "$output_file" "$log_file" &
    else
        (
            set -euo pipefail
            trap 'rm -f "$output_file"' ERR

            if [[ ! -s "$input_file" ]]; then
                echo "Missing or empty HepMC input: $input_file" >&2
                exit 1
            fi

            "$delphes_bin" "$rendered_card" "$output_file" "$input_file"

            python validate_delphes_output.py --hepmc "$input_file" --root "$output_file"
        ) > "$log_file" 2>&1 &
    fi
    pids+=("$!")
    logs+=("$log_file")

    job=$((job+1))
    if [[ -z "$stream_mg_tag" && "$job" == "$max_cpu_cores" ]]; then
        echo -e "\tStopping job submissions! Please wait for batch $batch to finish..."
        wait_workers
        job=0
        batch=$((batch+1))
    fi
done

wait_workers

echo -e "\tDelphes Simulation Done!"
