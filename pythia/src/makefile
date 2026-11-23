SUBMODULES_PATH=../../submodules/

CPP=g++
CXXFLAGS=-O2 -std=c++17 -pedantic -W -Wall -Wshadow -fPIC -pthread -w -DHEPMC2
PYTHIA_INC=-I$(SUBMODULES_PATH)/pythia-v8.312/include
PYTHIA_LIB=-L$(SUBMODULES_PATH)/pythia-v8.312/lib -Wl,-rpath,$(SUBMODULES_PATH)/pythia-v8.312/lib -lpythia8 -ldl
HEPMC_INC=-I$(SUBMODULES_PATH)/hepmc-v2.06.11/hepmc-install/include
HEPMC_LIB=-L$(SUBMODULES_PATH)/hepmc-v2.06.11/hepmc-install/lib -Wl,-rpath,$(SUBMODULES_PATH)/hepmc-v2.06.11/hepmc-install/lib -lHepMC

generate_dataset: generate_dataset.cpp
	@$(CPP) generate_dataset.cpp -o run_dataset $(CXXFLAGS) $(PYTHIA_INC) $(PYTHIA_LIB) $(HEPMC_INC) $(HEPMC_LIB)

clean:
	@rm -f run_dataset
