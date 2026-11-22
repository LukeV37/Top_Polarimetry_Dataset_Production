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
The training data will be written to a root file. 

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
