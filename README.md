## Quick Start
Clone the repo:
```
git clone --recursive git@github.com:LukeV37/Top_Polarimetry.git
```

Ensure dependencies are met, and install the submodules:
```
./build_submodules.sh
```

Please be patient while submodules build...

## How To Generate Datasets

### Madgraph
To run madgraph simulation, run the following
```
cd madgraph
./run.sh
./post_process.sh
```
The lhe file file is generated, and the training labels are written to root file. 

### Pythia
To shower in pythia, run the following
```
cd pythia
./run.sh
```
Pythia writes one HepMC file per run. Delphes reads those files and writes
reconstructed objects to ROOT.

### ROOT analysis selection (MVP)

Run the debug configuration from the repository root with ROOT (`root-config`)
and the built Delphes library available:

```bash
./run_job.sh debug.config
```

Set the bypass options in `debug.config` to skip earlier stages when processing
existing inputs. `run_job.sh` passes the analysis settings from the config to
the selection; a direct invocation of `analysis/run.sh` instead uses defaults
or explicitly supplied analysis environment variables.

This pairs Delphes and MadGraph label files by `_i` job suffix, applies the
selection, and writes `analysis/WS_<DELPHES_TAG>/dataset_selected_<DELPHES_TAG>_<i>.root`
and a per-job cutflow log under `logs/`. Each output has one `selected` TTree.
Set `ANALYSIS_OUTPUT_DIR` to an absolute path to write elsewhere. Selection
thresholds and array capacities are listed in `job.config` and `debug.config`;
the overlap threshold defaults to 1.5 and is independent of the jet radii
recorded in the rendered Delphes card. See the
[selection and branch reference](analysis/README.md) for
the cutflow, output fields, and overflow behavior. Files are not merged by
the runner; only merge successful outputs with identical settings and schemas.

## Build dependencies

On Ubuntu, install the system packages needed to configure and build the
submodules before running `./build_submodules.sh`:

```bash
sudo apt update
sudo apt install autoconf automake libtool
sudo apt install build-essential pkg-config
```

`autoconf`, `automake`, and `libtool` provide the GNU build/configuration
utilities used by the source builds. `build-essential` provides the compiler
and standard build tools (including `g++` and `make`); `pkg-config` helps build
scripts find installed libraries.

The build and dataset scripts also require:

- `python3` (and `python3-venv` if using `setup.sh` to create its virtual
  environment);
- `gfortran`;
- `gzip`;
- ROOT 6, installed separately and available through `root-config` (the apt
  commands above do not install ROOT).

The Pythia configuration summary should show the repository-built FastJet and
HepMC libraries, and ROOT if it is available in the environment. FastJet and
HepMC are built as submodules by `build_submodules.sh`.
