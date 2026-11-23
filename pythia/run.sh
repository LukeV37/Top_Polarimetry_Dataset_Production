#!/bin/bash

if [ -z "$4" ]; then
    echo "Must enter 4 arguments"
    echo -e "\t1: Dataset Tag (from MadGraph)"
    echo -e "\t2: Dataset Tag (for Pythia/HepMC output)"
    echo -e "\t3: Number of runs (from MadGraph)"
    echo -e "\t4: Max num cpu cores"
    exit 1
fi

MG_tag=$1
PY_tag=$2
runs=$3
max_cpu_cores=$4

mkdir -p "WS_${PY_tag}"
mkdir -p "WS_${PY_tag}/logs"

cd src || exit 1
make generate_dataset || exit 1

echo "Please be patient while Pythia performs hadronization..."

job=0
batch=1
pids=()
for (( i=0 ; i<$runs ; i++ ));
do
    echo -e "\t\tSubmitting job to shower run $i"
    ./run_dataset $MG_tag $PY_tag $i > "../WS_${PY_tag}/logs/dataset_${PY_tag}_${i}.log" 2>&1 &
    pids+=("$!")
    job=$((job+1))
    if [ $job == $max_cpu_cores ]; then
        echo -e "\tStopping jobs submissions! Please wait for batch $batch to finish..."
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
make clean || exit 1
cd .. || exit 1

echo -e "\tPythia Showering Done!"
