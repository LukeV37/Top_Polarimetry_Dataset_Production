#!/bin/bash
cd submodules/install_scripts
./clean_submodules.sh
./install_madgraph.sh
./install_hepmc.sh
./install_fastjet.sh
./install_pythia.sh
./install_delphes.sh
