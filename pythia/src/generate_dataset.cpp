#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        std::cout << "Error! Must enter 3 arguments and optionally an output path" << std::endl;
        std::cout << "1: Dataset Tag (from MadGraph)" << std::endl;
        std::cout << "2: Dataset Tag (for Pythia/HepMC output)" << std::endl;
        std::cout << "3: Run number (from MadGraph)" << std::endl;
        std::cout << "4: HepMC output path override (optional)" << std::endl;
        return 1;
    }

    const std::string in_tag = argv[1];
    const std::string out_tag = argv[2];
    const std::string run_num = argv[3];

    const std::string inputFile = "../../madgraph/pp_tt_semi_full_" + in_tag
        + "/Events/run_01_" + run_num + "/unweighted_events.lhe.gz";
    const std::string outputFile = argc == 5 ? argv[4] : "../WS_" + out_tag + "/dataset_" + out_tag
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

    std::ofstream output(outputFile);
    if (!output) {
        std::cerr << "Failed to open HepMC output file " << outputFile << std::endl;
        return 1;
    }
    HepMC::Pythia8ToHepMC toHepMC;
    auto writer = std::make_unique<HepMC::IO_GenEvent>(output);

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
            HepMC::GenEvent event;
            if (!toHepMC.fill_next_event(pythia, event)) {
                std::cerr << "HepMC conversion failed for run " << run_num
                          << " at event " << eventsWritten << std::endl;
                return 1;
            }
            writer->write_event(&event);
        } catch (const std::exception& ex) {
            std::cerr << "Exception while writing HepMC event for run " << run_num
                      << ": " << ex.what() << std::endl;
            return 1;
        }

        if (!output) {
            std::cerr << "HepMC output stream entered a failure state for run "
                      << run_num << " at event " << eventsWritten << std::endl;
            return 1;
        }

        ++eventsWritten;
    }

    writer.reset(); // Write the HepMC footer before closing and checking the stream.
    output.close();
    if (!output) {
        std::cerr << "Failed to finish HepMC output for run " << run_num << std::endl;
        return 1;
    }

    std::cout << "Successfully wrote " << eventsWritten << " events to "
              << outputFile << std::endl;
    if (argc == 5) {
        // Completion marker consumed by delphes/validate_delphes_output.py.
        std::cout << "HEPMC_EVENTS_WRITTEN=" << eventsWritten << std::endl;
    }
    return 0;
}
