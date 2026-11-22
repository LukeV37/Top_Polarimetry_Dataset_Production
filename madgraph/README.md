# MadGraph generation and truth labels

MadGraph is installed by the repository setup at
`submodules/mg5amcnlo-v3.5.17`. From the repository root, run it through
`./run_job.sh <config>`, or call `madgraph/run.sh` with the arguments shown by
the script.

## Directory layout

- `configs/` contains the process card, multi-run template, reference run card,
  and the custom hadronic-top pT cut.
- `include/` contains `LHEReader.py` and `TLorentz_Labels.py`. The reader accepts
  MadGraph's `.lhe.gz` files directly; it does not create an uncompressed LHE
  copy.
- MadGraph-generated process directories and event/ROOT files are runtime
  outputs, not source configuration.

## Hadronic-top pT cut

The configured `top_pT_cut` value is in GeV. The cut reconstructs the
parton-level hadronic top from external particles 3, 4, and 5 (`b j j`) and
rejects events when the transverse momentum of their summed four-momentum is
below the threshold. The particle ordering follows the current process card.

`configs/run_card.dat` carries the `min_hadronic_top_pt` entry, and
`madgraph/run.sh` sets both `custom_fcts` and `min_hadronic_top_pt` through the
existing `multi_run.tmp` command file. MadGraph then installs
`configs/hadronic_top_pt_cut.f` as the generated `dummy_cuts(P)` function. The
workflow does not patch generated `SubProcesses/cuts.f`.

The production and debug configs currently use `top_pT_cut=200`. The value is
passed through `run_job.sh` to `madgraph/run.sh`; changing it does not require
editing Fortran source.

`generation` can be `inclusive`, `first`, or `second`. The default inclusive
process keeps `t > b j j`; `first` uses `t > b u d~`; `second` uses
`t > b c s~`.
