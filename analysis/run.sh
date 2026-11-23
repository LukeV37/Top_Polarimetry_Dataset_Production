#!/bin/bash

set -euo pipefail

if [[ $# -ne 4 || ! $3 =~ ^[1-9][0-9]*$ || ! $4 =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: $0 MG_TAG DELPHES_TAG NUM_RUNS MAX_CPU_CORES" >&2
    exit 1
fi

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
MG_tag=$1
DELPHES_tag=$2
runs=$3
max_cpu_cores=$4
delphes_dir="../delphes/WS_${DELPHES_tag}"
labels_dir="../madgraph/pp_tt_semi_full_${MG_tag}"
workspace="${ANALYSIS_OUTPUT_DIR:-WS_${DELPHES_tag}}"
card="$delphes_dir/config/delphes_card_${DELPHES_tag}.tcl"

# Read the literal radius assignments for the log, not the dataset tag or the
# current config. The overlap threshold is configured independently.
read_radius() {
    awk -v key="$1" '
        $1 == "set" && $2 == key {count++; radius = $3; fields = NF}
        END {
            if (count != 1 || fields != 3 || radius !~ /^([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/ || radius <= 0) exit 1
            print radius
        }
    ' "$card" || { echo "Expected one positive $1 assignment in $card" >&2; return 1; }
}

large_r=$(read_radius LargeRJetRadius)
small_r=$(read_radius SmallRJetRadius)
for ((i = 0; i < runs; ++i)); do
    for input in "$delphes_dir/dataset_${DELPHES_tag}_${i}.root" "$labels_dir/labels_${MG_tag}_${i}.root"; do
        if [[ ! -s "$input" ]]; then
            echo "Missing or empty ROOT input: $input" >&2
            exit 1
        fi
    done
done

make run_selections
export ROOT_INCLUDE_PATH="$PWD/../submodules/delphes-v3.5.0/external${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
mkdir -p -- "$workspace/logs"

# Wait for every worker in the current batch; don't submit another batch after
# a failure. Each worker only publishes its ROOT file once it has succeeded.
pids=()
logs=()
wait_batch() {
    local i failed=0
    for i in "${!pids[@]}"; do
        if ! wait "${pids[i]}"; then
            echo "Selection job failed; see ${logs[i]}" >&2
            failed=1
        fi
    done
    pids=()
    logs=()
    return "$failed"
}

for ((i = 0; i < runs; ++i)); do
    stem="dataset_selected_${DELPHES_tag}_${i}"
    log="$workspace/logs/${stem}.log"
    echo "Submitting analysis for run $i ($log)"
    ./run_selections "$delphes_dir/dataset_${DELPHES_tag}_${i}.root" \
        "$labels_dir/labels_${MG_tag}_${i}.root" "$workspace/${stem}.root" \
        "$large_r" "$small_r" > "$log" 2>&1 &
    pids+=("$!")
    logs+=("$log")
    if ((${#pids[@]} == max_cpu_cores)); then
        wait_batch
    fi
done
wait_batch
echo "Analysis done: $runs selected ROOT files in $workspace"
