CXX ?= g++
ROOT_CONFIG ?= root-config
DELPHES_DIR ?= ../submodules/delphes-v3.5.0
DELPHES_PATH := $(abspath $(DELPHES_DIR))
CXXFLAGS ?= -O2 -Wall -Wextra

run_selections: selections.cxx $(DELPHES_PATH)/libDelphes.so
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $$($(ROOT_CONFIG) --cflags) \
	    -I$(DELPHES_PATH) -I$(DELPHES_PATH)/external $< -o $@ $(LDFLAGS) \
	    -L$(DELPHES_PATH) -Wl,-rpath,$(DELPHES_PATH) -lDelphes $$($(ROOT_CONFIG) --libs)

.PHONY: clean
clean:
	rm -f run_selections
