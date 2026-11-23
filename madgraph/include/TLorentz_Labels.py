import ROOT
from ROOT import TLorentzVector, TVector3
import numpy as np
import uproot
import awkward as ak
import sys

class Event:
    def __init__(self):
        # Inputs
        self.b = self.particle()
        self.u = self.particle()
        self.d = self.particle()
        self.bbar = self.particle()
        self.lep = self.particle()
        self.nu = self.particle()
        
        # Explicit, unnormalized four-vectors for the public label tree.
        self.top_lab = self.particle()
        self.top_ttCM = self.particle()
        self.top_topRest = self.particle()
        self.down_lab = self.particle()
        self.down_ttCM = self.particle()
        self.down_topRest = self.particle()
        self.b_lab = self.particle()
        self.b_ttCM = self.particle()
        self.b_topRest = self.particle()
        self.W_lab = self.particle()
        self.W_ttCM = self.particle()
        self.W_topRest = self.particle()
        self.costheta = []
        
    def __len__(self):
        length = len(self.b.px)
        if (len(self.u.px)!=length 
            or len(self.d.px)!=length 
            or len(self.bbar.px)!=length 
            or len(self.lep.px)!=length 
            or len(self.nu.px)!=length):
            raise Exception("Event does not contain uniquely identified quarks!"
                            +"Please refine selections.")
        return length
        
    class particle:
        def __init__(self):
            self.px = []
            self.py = []
            self.pz = []
            self.E = []
            self.pT = []
            self.eta = []
            self.phi = []
            self.deltaR = []
            self.deltaEta = []
            self.deltaPhi = []
            
        def add_feats(self, px, py, pz, E):
            self.px = np.concat((self.px, ak.to_numpy(px)))
            self.py = np.concat((self.py, ak.to_numpy(py)))
            self.pz = np.concat((self.pz, ak.to_numpy(pz)))
            self.E = np.concat((self.E, ak.to_numpy(E)))

        def initialize_four_vector(self, num_events):
            self.px = np.full(num_events, -999.0, dtype=np.float32)
            self.py = np.full(num_events, -999.0, dtype=np.float32)
            self.pz = np.full(num_events, -999.0, dtype=np.float32)
            self.E = np.full(num_events, -999.0, dtype=np.float32)

        def set_four_vector(self, event, vector):
            self.px[event] = vector.Px()
            self.py[event] = vector.Py()
            self.pz[event] = vector.Pz()
            self.E[event] = vector.E()
            
    def calc_labels(self):
        self.costheta_down = np.ones_like(self.b.px)*-999
        self.costheta_bottom = np.ones_like(self.b.px)*-999

        explicit_vectors = (
            self.top_lab, self.top_ttCM, self.top_topRest,
            self.down_lab, self.down_ttCM, self.down_topRest,
            self.b_lab, self.b_ttCM, self.b_topRest,
            self.W_lab, self.W_ttCM, self.W_topRest,
        )
        for vector in explicit_vectors:
            vector.initialize_four_vector(len(self.b.px))
        
        # Loop over events and calculate labels
        for event in range(len(self)):
            #Define particle's 4-vectors (px,py,pz,E)
            p_b = TLorentzVector(self.b.px[event],self.b.py[event],self.b.pz[event],self.b.E[event])
            p_d = TLorentzVector(self.d.px[event],self.d.py[event],self.d.pz[event],self.d.E[event])
            p_u = TLorentzVector(self.u.px[event],self.u.py[event],self.u.pz[event],self.u.E[event])

            p_bbar = TLorentzVector(self.bbar.px[event],self.bbar.py[event],self.bbar.pz[event],self.bbar.E[event])
            p_lm = TLorentzVector(self.lep.px[event],self.lep.py[event],self.lep.pz[event],self.lep.E[event])
            p_vl = TLorentzVector(self.nu.px[event],self.nu.py[event],self.nu.pz[event],self.nu.E[event])

            #Construct tops in lab frame by adding the corresponding decay products
            p_t = p_b + p_d + p_u 
            p_tbar = p_bbar + p_lm + p_vl
            # Reconstruct the hadronic W from its two quark daughters. The W
            # may not be present as a separate particle in the LHE record.
            p_W = p_u + p_d

            # Save actual lab-frame four-vectors before applying any boosts.
            self.top_lab.set_four_vector(event, p_t)
            self.b_lab.set_four_vector(event, p_b)
            self.down_lab.set_four_vector(event, p_d)
            self.W_lab.set_four_vector(event, p_W)

            #Construct Lorentz boost to t-tbar CM frame 
            to_ttbar_rest = -(p_t + p_tbar).BoostVector()

            #Boost vectors to t-tbar CM frame
            p_t.Boost(to_ttbar_rest)
            p_tbar.Boost(to_ttbar_rest)
            p_d.Boost(to_ttbar_rest)
            p_b.Boost(to_ttbar_rest)
            p_W.Boost(to_ttbar_rest)

            # Save actual ttbar-CM-frame four-vectors.
            self.top_ttCM.set_four_vector(event, p_t)
            self.down_ttCM.set_four_vector(event, p_d)
            self.b_ttCM.set_four_vector(event, p_b)
            self.W_ttCM.set_four_vector(event, p_W)

            #Top quark direction in t-tbar CM frame
            k_vect = p_t.Vect().Unit()

            #Define boost from t-tbar CM frame to top rest frame
            to_t_rest = -p_t.BoostVector()

            # Save the top four-vector in its own rest frame as well.
            p_t_topRest = TLorentzVector(p_t)
            p_t_topRest.Boost(to_t_rest)
            self.top_topRest.set_four_vector(event, p_t_topRest)

            #Boost the hadronic decay products to the top quark rest frame.
            p_d.Boost(to_t_rest)
            p_b.Boost(to_t_rest)
            p_W.Boost(to_t_rest)

            # Save actual, unnormalized top-rest-frame four-vectors.
            self.down_topRest.set_four_vector(event, p_d)
            self.b_topRest.set_four_vector(event, p_b)
            self.W_topRest.set_four_vector(event, p_W)

            #down quark direction in top rest frame
            d_vect = p_d.Vect().Unit()
            b_vect = p_b.Vect().Unit()

            #\cos\theta_d = k_vect\dot d_vect
            self.costheta_down[event] = k_vect.Dot(d_vect)
            self.costheta_bottom[event] = k_vect.Dot(b_vect)
            
    def write_ntuple(self, tag, num):
        with uproot.recreate("pp_tt_semi_full_"+tag+"/labels_"+tag+"_"+num+".root") as f:
            branches = {"costheta_down": "float32",
                        "costheta_bottom": "float32"}

            explicit_vectors = {
                "top": {"lab": self.top_lab,
                        "ttCM": self.top_ttCM,
                        "topRest": self.top_topRest},
                "bottom": {"lab": self.b_lab,
                           "ttCM": self.b_ttCM,
                           "topRest": self.b_topRest},
                "down": {"lab": self.down_lab,
                         "ttCM": self.down_ttCM,
                         "topRest": self.down_topRest},
                "W": {"lab": self.W_lab,
                      "ttCM": self.W_ttCM,
                      "topRest": self.W_topRest},
            }
            components = ("px", "py", "pz", "E")
            for object_name, frames in explicit_vectors.items():
                for frame_name in frames:
                    for component in components:
                        branches[f"{object_name}_{component}_{frame_name}"] = "float32"

            tree = f.mktree("labels", branches)

            data = {"costheta_down": self.costheta_down,
                    "costheta_bottom": self.costheta_bottom}

            for object_name, frames in explicit_vectors.items():
                for frame_name, vector in frames.items():
                    for component in components:
                        data[f"{object_name}_{component}_{frame_name}"] = getattr(vector, component)

            tree.extend(data)

# Get dataset tag
dataset_tag=str(sys.argv[1])
run_num=str(sys.argv[2])

# Read input ntuple
#print("\tReading lhe file...")
with uproot.open('pp_tt_semi_full_'+dataset_tag+'/hard_process_'+dataset_tag+'_'+run_num+'.root:events') as f:
    #print(f.keys())
    px = f['px'].array()
    py = f['py'].array()
    pz = f['pz'].array()
    E  = f['energy'].array()
    pid= f['pid'].array()
    m1 = f['mother1'].array()
    m2 = f['mother2'].array()
    status = f['status'].array()

# Initialize Event object
event = Event()

def add_one_particle_per_event(mask, particle, particle_name):
    """Keep the existing PID-mask approach and ensure event order is intact."""
    counts = np.asarray(ak.to_numpy(ak.sum(mask, axis=1)))
    bad_events = np.flatnonzero(counts != 1)
    if len(bad_events):
        examples = bad_events[:10].tolist()
        raise ValueError(
            f"Expected exactly one final-state {particle_name} per event; "
            f"found counts {counts[bad_events[:10]].tolist()} in entries {examples}"
        )

    particle.add_feats(
        np.ravel(px[mask]), np.ravel(py[mask]), np.ravel(pz[mask]), np.ravel(E[mask])
    )


# This is a fixed MadGraph process: identify its decay products by PID and
# final-state status, without adding ancestry-based selection.
final_state = status > 0
abs_pid = abs(pid)
particle_masks = (
    ((abs_pid == 1) | (abs_pid == 3)) & final_state,
    ((abs_pid == 2) | (abs_pid == 4)) & final_state,
    (pid == 5) & final_state,
    (pid == -5) & final_state,
    ((abs_pid == 11) | (abs_pid == 13)) & final_state,
    ((abs_pid == 12) | (abs_pid == 14)) & final_state,
)
particle_targets = (
    (event.d, "down-type quark"),
    (event.u, "up-type quark"),
    (event.b, "bottom quark"),
    (event.bbar, "antibottom quark"),
    (event.lep, "charged lepton"),
    (event.nu, "neutrino"),
)
for mask, (particle, particle_name) in zip(particle_masks, particle_targets):
    add_one_particle_per_event(mask, particle, particle_name)

# Calculate and write labels
#print("\tCalculating Labels...")
event.calc_labels()
#print("\tWriting ntuple...")
event.write_ntuple(dataset_tag,run_num)
#print("\tDone!")
