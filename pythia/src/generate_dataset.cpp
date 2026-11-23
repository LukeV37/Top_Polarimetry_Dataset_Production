#include <exception>
#include <ios>
#include <iostream>
#include <string>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Error! Must enter 3 arguments" << std::endl;
        std::cout << "1: Dataset Tag (from MadGraph)" << std::endl;
        std::cout << "2: Dataset Tag (for Pythia/HepMC output)" << std::endl;
        std::cout << "3: Run number (from MadGraph)" << std::endl;
        return 1;
    }

    const std::string in_tag = argv[1];
    const std::string out_tag = argv[2];
    const std::string run_num = argv[3];

    const std::string inputFile = "../../madgraph/pp_tt_semi_full_" + in_tag
        + "/Events/run_01_" + run_num + "/unweighted_events.lhe.gz";
    const std::string outputFile = "../WS_" + out_tag + "/dataset_" + out_tag
        + "_" + run_num + ".hepmc";

    std::cout << "Input LHE: " << inputFile << std::endl;
    std::cout << "Output HepMC: " << outputFile << std::endl;

    Pythia8::Pythia pythia;
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = " + inputFile);
    pythia.readString("Next:numberCount = 1000");

    if (!pythia.init()) {
        std::cerr << "Pythia failed to initialize for input " << inputFile << std::endl;
        return 1;
    }

    Pythia8::Pythia8ToHepMC toHepMC;
    if (!toHepMC.setNewFile(outputFile)) {
        std::cerr << "Failed to open HepMC output file " << outputFile << std::endl;
        return 1;
    }

    int eventsWritten = 0;

    while (true) {
        if (!pythia.next()) {
            if (pythia.info.atEndOfFile()) {
                break;
            }

            std::cerr << "Pythia failed before EOF for run " << run_num
                      << " after " << eventsWritten << " successful events." << std::endl;
            return 1;
        }

        try {
            if (!toHepMC.writeNextEvent(pythia)) {
                std::cerr << "HepMC conversion failed for run " << run_num
                          << " at event " << eventsWritten << std::endl;
                return 1;
            }
        } catch (const std::exception& ex) {
            std::cerr << "Exception while writing HepMC event for run " << run_num
                      << ": " << ex.what() << std::endl;
            return 1;
        }

        const int streamState = toHepMC.output().rdstate();
        if ((streamState & std::ios::failbit) || (streamState & std::ios::badbit)) {
            std::cerr << "HepMC output stream entered a failure state for run "
                      << run_num << " at event " << eventsWritten << std::endl;
            return 1;
        }

        ++eventsWritten;
    }

    std::cout << "Successfully wrote " << eventsWritten << " events to "
              << outputFile << std::endl;
    return 0;
}
