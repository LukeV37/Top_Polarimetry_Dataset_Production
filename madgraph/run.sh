#!/bin/bash

if [[ $# -ne 8 ]]; then
    echo "Must enter 8 arguments"
    echo -e "\t1: Dataset Tag e.g. unpolarized_10k"
    echo -e "\t2: Polarization: L, R, U"
    echo -e "\t3: Generation: inclusive, first, second"
    echo -e "\t4: Number of Runs"
    echo -e "\t5: Number of Events per Run"
    echo -e "\t6: Max num cpu cores"
    echo -e "\t7: Initial seed"
    echo -e "\t8: Hadronic-top pT threshold in GeV"
    exit 1
fi

dataset_tag=$1
polarization=$2
generation=$3
num_runs=$4
num_events_per_run=$5
max_cpu_cores=$6
seed=$7
top_pT_cut=$8

# Error handling before launching event generation
set -euo pipefail
if [[ ! $top_pT_cut =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v cut="$top_pT_cut" 'BEGIN { exit !(cut > 0) }'; then
    echo "Hadronic-top pT threshold must be a positive number in GeV: $top_pT_cut" >&2
    exit 1
fi
case "$generation" in
    inclusive|first|second) ;;
    *)
        echo "Generation must be one of: inclusive, first, second. Got: $generation" >&2
        exit 1
        ;;
esac
python -c "import ROOT; import uproot; import awkward"

custom_cut_file=$(realpath configs/hadronic_top_pt_cut.f)
if [[ ! -f "$custom_cut_file" ]]; then
    echo "Custom cut file not found: $custom_cut_file" >&2
    exit 1
fi

###########################
### MadGraph Generation ###
###########################

# Edit process card
sed "s/output pp_tt_semi_full/output pp_tt_semi_full_${dataset_tag}/" configs/proc_card_mg5.dat > proc_card.tmp
if [ "$polarization" = "L" ]; then
  sed -i "s/generate p p > t t~, t > b j j, t~ > b~ l- vl~/generate p p > t{L} t~, t > b j j, t~ > b~ l- vl~/" proc_card.tmp
fi
if [ "$polarization" = "R" ]; then
  sed -i "s/generate p p > t t~, t > b j j, t~ > b~ l- vl~/generate p p > t{R} t~, t > b j j, t~ > b~ l- vl~/" proc_card.tmp
fi
if [ "$polarization" = "down" ]; then
  sed -i "s/generate p p > t t~, t > b j j, t~ > b~ l- vl~/generate p p > t t~, t > b u d~, t~ > b~ l- vl~/" proc_card.tmp
fi
if [ "$polarization" = "charm" ]; then
  sed -i "s/generate p p > t t~, t > b j j, t~ > b~ l- vl~/generate p p > t t~, t > b c s~, t~ > b~ l- vl~/" proc_card.tmp
fi
if [ "$generation" = "first" ]; then
  sed -i "s/t > b j j/t > b u d~/" proc_card.tmp
elif [ "$generation" = "second" ]; then
  sed -i "s/t > b j j/t > b c s~/" proc_card.tmp
fi

# Prepare the multi_run command with the requested run count, event count, and seed.
sed "s/multi_run.*/multi_run $num_runs/" configs/multi_run.config > multi_run.tmp
sed -i "s/set nevents.*/set nevents $num_events_per_run/" multi_run.tmp
sed -i "s/set iseed.*/set iseed $seed/" multi_run.tmp
printf "  set custom_fcts ['%s']\n" "$custom_cut_file" >> multi_run.tmp
printf "  set min_hadronic_top_pt %s\n" "$top_pT_cut" >> multi_run.tmp

# Run MadGraph 3.5.17 from the tarball install
../submodules/mg5amcnlo-v3.5.17/bin/mg5_aMC proc_card.tmp

# Install the checked-in run card template before event generation.
process_dir="./pp_tt_semi_full_${dataset_tag}"
cp configs/run_card.dat "$process_dir/Cards/run_card.dat"
grep -i "min_hadronic_top_pt" "$process_dir/Cards/run_card.dat" >/dev/null

# Generate the events.
echo "Please be patient while MadGraph generates processes..."
"$process_dir/bin/madevent" multi_run.tmp | tee "MadGraph_${dataset_tag}.log"

# Verify MadGraph applied the runtime run-card overrides and installed the custom
# dummy_cuts hook instead of patching cuts.f.
banner_file="$process_dir/Events/run_01_banner.txt"
grep -F "$custom_cut_file" "$banner_file" >/dev/null
grep -i "min_hadronic_top_pt" "$banner_file" >/dev/null
grep -i "min_hadronic_top_pt" "$process_dir/SubProcesses/dummy_fct.f" >/dev/null
if grep -q "CUSTOM_HADRONIC_TOP_PT_CUT" "$process_dir/SubProcesses/cuts.f"; then
    echo "Unexpected generated cuts.f patch marker found; custom_fcts path failed" >&2
    exit 1
fi
for (( run=0 ; run<num_runs ; run++ )); do
    lhe_file="$process_dir/Events/run_01_$run/unweighted_events.lhe.gz"
    if [[ ! -s "$lhe_file" ]]; then
        echo "MadGraph did not produce expected LHE file: $lhe_file" >&2
        exit 1
    fi
done

# Clean up workspace (generated automatically by madgraph binary)
rm -f py.py
rm -f MG5_debug
rm -f ME5_debug
rm -f -- *.tmp

#######################
### Post Processing ###
#######################

function calc_labels {
    tag=$1
    run=$2

    # Read the compressed LHE directly; LHEReader streams it through gzip.
    #echo -e "\tConverting lhe file to root format..."
    python include/LHEReader.py --input "pp_tt_semi_full_${tag}/Events/run_01_$run/unweighted_events.lhe.gz" --output "pp_tt_semi_full_${tag}/hard_process_${tag}_$run.root"

    # Calculate labels
    python include/TLorentz_Labels.py "$tag" "$run"
}

job=0
batch=1
for (( i=0 ; i<$num_runs ; i++ ));
do
    echo -e "\t\tSubmitting job to calc labels: $i"
    calc_labels "$dataset_tag" "$i" &
    job=$((job+1))
    if [ $job == $max_cpu_cores ]; then
        echo -e "\tStopping jobs submissions! Please wait for batch $batch to finish..."
        wait || exit 1
        job=0
        batch=$((batch+1))
    fi
done
wait || exit 1

for (( run=0 ; run<num_runs ; run++ )); do
    for stage in hard_process labels; do
        output_file="pp_tt_semi_full_${dataset_tag}/${stage}_${dataset_tag}_${run}.root"
        if [[ ! -s "$output_file" ]]; then
            echo "Missing or empty ROOT output: $output_file" >&2
            exit 1
        fi
    done
done

echo
echo -e "\tMadGraph Generation Done!"
