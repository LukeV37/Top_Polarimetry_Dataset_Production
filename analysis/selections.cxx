#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <unistd.h>

#include "TBranchElement.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TTree.h"
#include "classes/DelphesClasses.h"

namespace {
constexpr int padding = -999;
constexpr double twoPi = 6.2831853071795864769;

double number(const char *name, const char *text)
{
    char *end = nullptr;
    const double value = std::strtod(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(value) || value < 0)
        throw std::runtime_error(std::string(name) + " must be a finite nonnegative number");
    return value;
}

double setting(const char *name, const char *fallback)
{
    const char *text = std::getenv(name);
    const double value = number(name, text ? text : fallback);
    std::cout << name << '=' << value << '\n';
    return value;
}

int capacity(const char *name, const char *fallback)
{
    const double value = setting(name, fallback);
    if (value < 1 || value > 100000 || std::floor(value) != value)
        throw std::runtime_error(std::string(name) + " must be an integer in [1, 100000]");
    return static_cast<int>(value);
}

struct Config {
    double electronPt = setting("ELECTRON_PT_MIN", "30");
    double electronEta = setting("ELECTRON_ABS_ETA_MAX", "3");
    double muonPt = setting("MUON_PT_MIN", "30");
    double muonEta = setting("MUON_ABS_ETA_MAX", "3");
    double met = setting("MET_MIN", "30");
    double largeJetPt = setting("LARGE_R_JET_PT_MIN", "250");
    double largeJetEta = setting("LARGE_R_JET_ABS_ETA_MAX", "3");
    double overlapDR = setting("JET_OVERLAP_DR", "1.5");
    double bJetPt = setting("LEPTONIC_B_JET_PT_MIN", "30");
    double bJetEta = setting("LEPTONIC_B_JET_ABS_ETA_MAX", "3");
    double bJetDRMin = setting("LEPTON_B_JET_DR_MIN", "0.3");
    double bJetDRMax = setting("LEPTON_B_JET_DR_MAX", "3.0");
    int maxConstituents = capacity("MAX_LARGE_R_CONSTITUENTS", "100");
    int maxSmallJets = capacity("MAX_SMALL_R_JETS", "10");
};

TTree &inputTree(TFile &file, const char *name)
{
    TTree *tree = nullptr;
    if (file.IsOpen() && !file.IsZombie() && !file.TestBit(TFile::kRecovered))
        file.GetObject(name, tree);
    if (!tree)
        throw std::runtime_error(std::string("Cannot read tree ") + name + " in " + file.GetName());
    return *tree;
}

void readBranch(TBranch &branch, Long64_t entry)
{
    if (branch.GetEntry(entry) <= 0)
        throw std::runtime_error(std::string("Cannot read branch ") + branch.GetName());
}

// Delphes collections are split TClonesArrays. Bind them as in ExRootTreeReader,
// but check each branch read rather than only the entry index.
class Collection {
public:
    Collection(TTree &tree, const char *name, TClass *type) : objects_(type), address_(&objects_)
    {
        branch_ = dynamic_cast<TBranchElement *>(tree.GetBranch(name));
        if (!branch_ || !branch_->GetClonesName() || std::string(branch_->GetClonesName()) != type->GetName())
            throw std::runtime_error(std::string("Missing or incompatible Delphes collection: ") + name);
        branch_->SetAddress(&address_);
    }

    ~Collection() { branch_->ResetAddress(); }
    Collection(const Collection &) = delete;
    Collection &operator=(const Collection &) = delete;

    void read(Long64_t entry) { readBranch(*branch_, entry); }
    int size() const { return objects_.GetEntriesFast(); }
    const TObject *at(int index) const { return objects_.At(index); }

    template <typename T>
    const T &get(int index) const
    {
        const auto *object = static_cast<const T *>(at(index));
        if (!object) throw std::runtime_error(std::string("Null object in ") + branch_->GetName());
        return *object;
    }

private:
    TClonesArray objects_;
    TClonesArray *address_;
    TBranchElement *branch_;
};

struct Label {
    std::string name;
    Float_t value = 0;
    TBranch *branch = nullptr;
};

void outputBranch(TTree &tree, const std::string &name, void *address, char type, int slots = 0)
{
    const auto leaf = name + (slots ? "[" + std::to_string(slots) + "]" : "") + '/' + type;
    if (!tree.Branch(name.c_str(), address, leaf.c_str()))
        throw std::runtime_error("Cannot create output branch " + name);
}

std::vector<Label> labelsForMvp()
{
    std::vector<Label> labels{{"costheta_down"}, {"costheta_bottom"}};
    for (const auto *particle : {"top", "bottom", "down", "W"})
        for (const auto *frame : {"lab", "ttCM", "topRest"})
            for (const auto *component : {"px", "py", "pz", "E"})
                labels.push_back({std::string(particle) + '_' + component + '_' + frame});
    return labels;
}

void bindLabel(TTree &input, TTree &output, Label &label)
{
    auto *branch = input.GetBranch(label.name.c_str());
    auto *leaf = branch ? branch->GetLeaf(label.name.c_str()) : nullptr;
    if (!leaf || std::string(leaf->GetTypeName()) != "Float_t" || leaf->GetLeafCount() ||
        leaf->GetLenStatic() != 1 || branch->GetNleaves() != 1 ||
        input.SetBranchAddress(label.name.c_str(), &label.value, &label.branch) != TTree::kMatch)
        throw std::runtime_error("Missing or incompatible MadGraph label: " + label.name);
    outputBranch(output, label.name, &label.value, 'F');
}

struct Counts {
    Int_t total = 0, stored = 0, overflow = 0;

    void bind(TTree &tree, const std::string &name)
    {
        outputBranch(tree, "num_" + name, &total, 'I');
        outputBranch(tree, "num_" + name + "_stored", &stored, 'I');
        outputBranch(tree, "num_" + name + "_overflow", &overflow, 'I');
    }

    void update(int count, int slots)
    {
        total = count;
        stored = std::min(total, slots);
        overflow = total - stored;
    }
};

// The buffer size never changes after binding: the ROOT leaf is a fixed array.
template <typename T>
class FixedArray {
public:
    FixedArray(TTree &tree, const char *name, int slots, char type) : values_(slots, padding)
    {
        outputBranch(tree, name, values_.data(), type, slots);
    }

    void padFrom(int first) { std::fill(values_.begin() + first, values_.end(), padding); }
    T &operator[](int index) { return values_[index]; }

private:
    std::vector<T> values_;
};

struct Constituent {
    Float_t pt, eta, phi, energy;
    Int_t type, charge, pid;
    Float_t d0, dz, d0Error, dzError;
};

std::vector<Constituent> getConstituents(const Jet &jet, const Collection &tracks,
                                         const Collection &photons, const Collection &neutrals)
{
    // Reject references that do not resolve to an energy-flow object in this entry.
    std::unordered_set<const TObject *> current;
    for (const auto *collection : {&tracks, &photons, &neutrals})
        for (int i = 0; i < collection->size(); ++i)
            current.insert(collection->at(i));

    std::vector<Constituent> result;
    for (int i = 0; i < jet.Constituents.GetEntriesFast(); ++i) {
        const auto *object = jet.Constituents.At(i);
        if (!object || !current.count(object))
            throw std::runtime_error("Unresolved large-R constituent reference");
        if (const auto *track = dynamic_cast<const Track *>(object)) {
            result.push_back({track->PT, track->Eta, track->Phi, static_cast<Float_t>(track->P4().E()),
                              1, track->Charge, track->PID, track->D0, track->DZ,
                              track->ErrorD0, track->ErrorDZ});
        } else if (const auto *tower = dynamic_cast<const Tower *>(object)) {
            result.push_back({tower->ET, tower->Eta, tower->Phi, tower->E, 0, padding, padding,
                              padding, padding, padding, padding});
        } else {
            throw std::runtime_error("Unexpected large-R constituent class");
        }
        if (!std::isfinite(result.back().pt))
            throw std::runtime_error("Non-finite large-R constituent pT");
    }
    std::stable_sort(result.begin(), result.end(), [](const auto &a, const auto &b) { return a.pt > b.pt; });
    return result;
}

double deltaR(double eta1, double phi1, double eta2, double phi2)
{
    return std::hypot(eta1 - eta2, std::remainder(phi1 - phi2, twoPi));
}

bool tagged(const Jet &jet) { return (jet.BTag & 1u) != 0; }

struct Kinematics {
    Float_t pt = 0, eta = 0, phi = 0, mass = 0;

    void bind(TTree &tree, const std::string &name)
    {
        outputBranch(tree, name + "_pt", &pt, 'F');
        outputBranch(tree, name + "_eta", &eta, 'F');
        outputBranch(tree, name + "_phi", &phi, 'F');
        outputBranch(tree, name + "_mass", &mass, 'F');
    }
};

void select(TTree &delphes, TTree &truth, TFile &output, const Config &cfg)
{
    const Long64_t entries = delphes.GetEntries();
    if (entries != truth.GetEntries())
        throw std::runtime_error("Entry-count mismatch: Delphes=" + std::to_string(entries) +
                                 ", labels=" + std::to_string(truth.GetEntries()));

    Collection events(delphes, "Event", HepMCEvent::Class());
    Collection electrons(delphes, "Electron", Electron::Class());
    Collection muons(delphes, "Muon", Muon::Class());
    Collection missingET(delphes, "MissingET", MissingET::Class());
    Collection largeJets(delphes, "LargeRJet", Jet::Class());
    Collection smallJets(delphes, "SmallRJet", Jet::Class());
    Collection tracks(delphes, "EFlowTrack", Track::Class());
    Collection photons(delphes, "EFlowPhoton", Tower::Class());
    Collection neutrals(delphes, "EFlowNeutralHadron", Tower::Class());

    output.cd();
    auto labels = labelsForMvp();
    TTree selected("selected", "Selected Delphes events with MadGraph truth labels");
    const int nc = cfg.maxConstituents, nj = cfg.maxSmallJets;
    selected.GetUserInfo()->Add(new TParameter<int>("MAX_LARGE_R_CONSTITUENTS", nc));
    selected.GetUserInfo()->Add(new TParameter<int>("MAX_SMALL_R_JETS", nj));
    selected.GetUserInfo()->Add(new TParameter<int>("padding_value", padding));
    selected.GetUserInfo()->Add(new TParameter<double>("JET_OVERLAP_DR", cfg.overlapDR));
    for (auto &label : labels) bindLabel(truth, selected, label);

    Long64_t sourceEntry = 0;
    Float_t leptonPt = 0, leptonEta = 0, leptonPhi = 0, isolation = 0, isolationRho = 0;
    Float_t met = 0, metPhi = 0;
    Int_t leptonCharge = 0, leptonFlavor = 0;
    Kinematics largeJet;
    Counts constituentCounts, smallJetCounts;
    outputBranch(selected, "source_entry", &sourceEntry, 'L');
    outputBranch(selected, "lepton_pt", &leptonPt, 'F');
    outputBranch(selected, "lepton_eta", &leptonEta, 'F');
    outputBranch(selected, "lepton_phi", &leptonPhi, 'F');
    outputBranch(selected, "lepton_charge", &leptonCharge, 'I');
    outputBranch(selected, "lepton_flavor", &leptonFlavor, 'I');
    outputBranch(selected, "lepton_isolation", &isolation, 'F');
    outputBranch(selected, "lepton_isolation_rho_corr", &isolationRho, 'F');
    outputBranch(selected, "met", &met, 'F');
    outputBranch(selected, "met_phi", &metPhi, 'F');
    largeJet.bind(selected, "large_r_jet");
    constituentCounts.bind(selected, "constituents");
    smallJetCounts.bind(selected, "small_r_jets");

    FixedArray<Float_t> cPt(selected, "constituent_pt", nc, 'F'), cEta(selected, "constituent_eta", nc, 'F'),
        cPhi(selected, "constituent_phi", nc, 'F'), cE(selected, "constituent_E", nc, 'F'),
        d0(selected, "constituent_d0", nc, 'F'), dz(selected, "constituent_dz", nc, 'F'),
        d0Error(selected, "constituent_d0_error", nc, 'F'), dzError(selected, "constituent_dz_error", nc, 'F');
    FixedArray<Int_t> cType(selected, "constituent_type", nc, 'I'), cCharge(selected, "constituent_charge", nc, 'I'),
        cPid(selected, "constituent_pid", nc, 'I');
    FixedArray<Float_t> jPt(selected, "small_r_jet_pt", nj, 'F'), jEta(selected, "small_r_jet_eta", nj, 'F'),
        jPhi(selected, "small_r_jet_phi", nj, 'F'), jMass(selected, "small_r_jet_mass", nj, 'F');
    FixedArray<Int_t> jTag(selected, "small_r_jet_btag", nj, 'I');

    Long64_t passLepton = 0, passMet = 0, passLargeJet = 0;
    Long64_t acceptedElectrons = 0, acceptedMuons = 0;
    Long64_t constituentOverflow = 0, smallJetOverflow = 0;
    for (sourceEntry = 0; sourceEntry < entries; ++sourceEntry) {
        try {
            events.read(sourceEntry);
            if (events.size() != 1 || events.get<HepMCEvent>(0).Number != sourceEntry)
                throw std::runtime_error("Non-sequential Delphes Event.Number");
            // Read both inputs even for entries that will be rejected by physics cuts.
            for (auto &label : labels) readBranch(*label.branch, sourceEntry);
            for (auto *collection : {&electrons, &muons, &missingET, &largeJets, &smallJets, &tracks, &photons, &neutrals})
                collection->read(sourceEntry);
            if (missingET.size() != 1) throw std::runtime_error("Expected one MissingET object");

            // Use the highest-pT qualifying lepton across both collections.
            // Strict '>' preserves collection order for ties, with electrons first.
            leptonFlavor = 0;
            auto consider = [&](const auto &candidate, int flavor, double minPt, double maxEta) {
                if (candidate.PT >= minPt && std::abs(candidate.Eta) <= maxEta &&
                    (!leptonFlavor || candidate.PT > leptonPt)) {
                    leptonPt = candidate.PT;
                    leptonEta = candidate.Eta;
                    leptonPhi = candidate.Phi;
                    leptonCharge = candidate.Charge;
                    leptonFlavor = flavor;
                    isolation = candidate.IsolationVar;
                    isolationRho = candidate.IsolationVarRhoCorr;
                }
            };
            for (int i = 0; i < electrons.size(); ++i) {
                consider(electrons.get<Electron>(i), 11, cfg.electronPt, cfg.electronEta);
            }
            for (int i = 0; i < muons.size(); ++i) {
                consider(muons.get<Muon>(i), 13, cfg.muonPt, cfg.muonEta);
            }
            if (!leptonFlavor) continue;
            ++passLepton;

            const auto &missing = missingET.get<MissingET>(0);
            if (!(missing.MET >= cfg.met)) continue;
            ++passMet;

            int numLarge = 0;
            const Jet *jet = nullptr;
            for (int i = 0; i < largeJets.size(); ++i) {
                const auto &candidate = largeJets.get<Jet>(i);
                if (candidate.PT >= cfg.largeJetPt) {
                    ++numLarge;
                    jet = &candidate;
                }
            }
            if (numLarge != 1 || !(std::abs(jet->Eta) <= cfg.largeJetEta)) continue;
            if (!std::isfinite(jet->Phi) || !std::isfinite(leptonPhi))
                throw std::runtime_error("Non-finite lepton or large-R jet azimuth");
            ++passLargeJet;
            auto constituents = getConstituents(*jet, tracks, photons, neutrals);

            bool hasB = false;
            std::vector<const Jet *> nonOverlapping;
            for (int i = 0; i < smallJets.size(); ++i) {
                const auto &small = smallJets.get<Jet>(i);
                if (!std::isfinite(small.PT) || !std::isfinite(small.Eta) || !std::isfinite(small.Phi))
                    throw std::runtime_error("Non-finite small-R jet kinematics");
                if (deltaR(jet->Eta, jet->Phi, small.Eta, small.Phi) < cfg.overlapDR) continue;
                nonOverlapping.push_back(&small);
                const double dr = deltaR(leptonEta, leptonPhi, small.Eta, small.Phi);
                if (tagged(small) && small.PT > cfg.bJetPt && std::abs(small.Eta) < cfg.bJetEta &&
                    dr > cfg.bJetDRMin && dr < cfg.bJetDRMax) hasB = true;
            }
            if (!hasB) continue;
            std::stable_sort(nonOverlapping.begin(), nonOverlapping.end(),
                             [](const Jet *a, const Jet *b) { return a->PT > b->PT; });

            met = missing.MET;
            metPhi = missing.Phi;
            largeJet = {jet->PT, jet->Eta, jet->Phi, jet->Mass};

            constituentCounts.update(static_cast<int>(constituents.size()), nc);
            constituentOverflow += constituentCounts.overflow;
            for (int i = 0; i < constituentCounts.stored; ++i) {
                const auto &c = constituents[i];
                cPt[i] = c.pt; cEta[i] = c.eta; cPhi[i] = c.phi; cE[i] = c.energy;
                cType[i] = c.type; cCharge[i] = c.charge; cPid[i] = c.pid;
                d0[i] = c.d0; dz[i] = c.dz; d0Error[i] = c.d0Error; dzError[i] = c.dzError;
            }
            for (auto *field : {&cPt, &cEta, &cPhi, &cE, &d0, &dz, &d0Error, &dzError})
                field->padFrom(constituentCounts.stored);
            for (auto *field : {&cType, &cCharge, &cPid}) field->padFrom(constituentCounts.stored);

            smallJetCounts.update(static_cast<int>(nonOverlapping.size()), nj);
            smallJetOverflow += smallJetCounts.overflow;
            for (int i = 0; i < smallJetCounts.stored; ++i) {
                const auto &small = *nonOverlapping[i];
                jPt[i] = small.PT; jEta[i] = small.Eta; jPhi[i] = small.Phi; jMass[i] = small.Mass;
                jTag[i] = tagged(small);
            }
            for (auto *field : {&jPt, &jEta, &jPhi, &jMass}) field->padFrom(smallJetCounts.stored);
            jTag.padFrom(smallJetCounts.stored);

            if (selected.Fill() < 0) throw std::runtime_error("Cannot fill selected tree");
            if (leptonFlavor == 11) ++acceptedElectrons;
            else ++acceptedMuons;
        } catch (const std::exception &error) {
            throw std::runtime_error("Entry " + std::to_string(sourceEntry) + ": " + error.what());
        }
    }

    if (selected.Write() <= 0 || output.TestBit(TFile::kWriteError))
        throw std::runtime_error("Cannot write selected tree");
    truth.ResetBranchAddresses();
    std::cout << "Input events: " << entries << '\n'
              << "Passed lepton: " << passLepton << '\n'
              << "Passed MET: " << passMet << '\n'
              << "Passed large-R jet: " << passLargeJet << '\n'
              << "Passed leptonic b jet (accepted): " << selected.GetEntries() << '\n'
              << "Accepted electrons: " << acceptedElectrons << '\n'
              << "Accepted muons: " << acceptedMuons << '\n'
              << "Constituents dropped by capacity: " << constituentOverflow << '\n'
              << "Small-R jets dropped by capacity: " << smallJetOverflow << '\n';
}
} // namespace

int main(int argc, char **argv)
{
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " DELPHES.root LABELS.root OUTPUT.root R_LARGE R_SMALL\n";
        return 1;
    }

    std::string temporary;
    try {
        const double largeR = number("R_LARGE", argv[4]);
        const double smallR = number("R_SMALL", argv[5]);
        if (smallR <= 0 || largeR <= smallR || !std::isfinite(largeR + smallR))
            throw std::runtime_error("Require finite jet radii with R_LARGE > R_SMALL > 0");
        std::cout << std::setprecision(16) << "Delphes: " << argv[1] << "\nLabels: " << argv[2]
                  << "\nR_LARGE=" << largeR << "\nR_SMALL=" << smallR << '\n';
        const Config cfg;
        if (cfg.bJetDRMin >= cfg.bJetDRMax)
            throw std::runtime_error("Require LEPTON_B_JET_DR_MIN < LEPTON_B_JET_DR_MAX");

        TFile delphesFile(argv[1], "READ"), truthFile(argv[2], "READ");
        auto &delphes = inputTree(delphesFile, "Delphes");
        auto &truth = inputTree(truthFile, "labels");
        if (std::filesystem::exists(argv[3]) && (std::filesystem::equivalent(argv[3], argv[1]) ||
                                               std::filesystem::equivalent(argv[3], argv[2])))
            throw std::runtime_error("Output would overwrite an input ROOT file");

        // Publish only a complete output; a failed rerun leaves any previous result intact.
        std::string pattern = std::string(argv[3]) + ".tmp.XXXXXX";
        const int fd = mkstemp(pattern.data());
        if (fd == -1) throw std::runtime_error("Cannot create temporary output next to " + std::string(argv[3]));
        temporary = std::move(pattern);
        close(fd);
        TFile output(temporary.c_str(), "RECREATE");
        if (!output.IsOpen() || output.IsZombie()) throw std::runtime_error("Cannot open output ROOT file");
        select(delphes, truth, output, cfg);
        output.Close();
        if (output.TestBit(TFile::kWriteError)) throw std::runtime_error("Error closing output file");
        std::filesystem::rename(temporary, argv[3]);
        std::cout << "Wrote " << argv[3] << '\n';
    } catch (const std::exception &error) {
        if (!temporary.empty()) {
            std::error_code ignored;
            std::filesystem::remove(temporary, ignored);
        }
        std::cerr << "Selection failed: " << error.what() << '\n';
        return 1;
    }
}
