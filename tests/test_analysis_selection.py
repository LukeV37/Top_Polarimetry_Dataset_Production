"""Validate selection against existing debug inputs without regenerating events."""

import math
import os
from pathlib import Path
import re
import shutil
import subprocess

import awkward as ak
import numpy as np
import pytest
import uproot


REPO = Path(__file__).resolve().parents[1]
ANALYSIS = REPO / "analysis"
DELPHES = REPO / "submodules/delphes-v3.5.0"
ANALYSIS_ENV_KEYS = {
    "ANALYSIS_OUTPUT_DIR", "ELECTRON_PT_MIN", "ELECTRON_ABS_ETA_MAX", "MUON_PT_MIN",
    "MUON_ABS_ETA_MAX", "MET_MIN", "LARGE_R_JET_PT_MIN", "LARGE_R_JET_ABS_ETA_MAX",
    "JET_OVERLAP_DR", "LEPTONIC_B_JET_PT_MIN", "LEPTONIC_B_JET_ABS_ETA_MAX",
    "LEPTON_B_JET_DR_MIN", "LEPTON_B_JET_DR_MAX", "MAX_LARGE_R_CONSTITUENTS", "MAX_SMALL_R_JETS",
}


@pytest.fixture(scope="module")
def root_build():
    if not shutil.which("root-config") or not (DELPHES / "libDelphes.so").is_file():
        pytest.skip("Requires ROOT and the built Delphes library")
    subprocess.run(["make", "run_selections"], cwd=ANALYSIS, check=True)


@pytest.fixture(scope="module")
def debug_inputs(root_build):
    mg, delphes, runs = subprocess.check_output(
        ["bash", "-c", 'source debug.config && printf "%s\n" "$MG_tag" "$DELPHES_tag" "$num_runs"'],
        cwd=REPO, text=True,
    ).splitlines()
    workspace = REPO / "delphes" / f"WS_{delphes}"
    reco = [workspace / f"dataset_{delphes}_{i}.root" for i in range(int(runs))]
    labels = [REPO / "madgraph" / f"pp_tt_semi_full_{mg}" / f"labels_{mg}_{i}.root" for i in range(int(runs))]
    card = workspace / "config" / f"delphes_card_{delphes}.tcl"
    if not all(path.is_file() for path in [*reco, *labels, card]):
        pytest.skip("Requires the existing debug.config Delphes outputs, card, and MadGraph labels")
    radii = {parts[1]: float(parts[2]) for line in card.read_text().splitlines()
             if len(parts := line.split()) == 3 and parts[0] == "set" and
             parts[1] in ("LargeRJetRadius", "SmallRJetRadius")}
    return mg, delphes, reco, labels, radii["LargeRJetRadius"], radii["SmallRJetRadius"]


def analysis_env(**settings):
    env = {key: val for key, val in os.environ.items() if key not in ANALYSIS_ENV_KEYS}
    env["ROOT_INCLUDE_PATH"] = f'{DELPHES / "external"}:{env.get("ROOT_INCLUDE_PATH", "")}'
    env.update(settings)
    return env


def records(tree, branch, *fields):
    fields = {field: tree[f"{branch}/{branch}.{field}"].array() for field in fields}
    return ak.to_list(ak.zip(fields, depth_limit=2))


def dr(a, b):
    dphi = (a["Phi"] - b["Phi"] + math.pi) % (2 * math.pi) - math.pi
    return math.hypot(a["Eta"] - b["Eta"], dphi)


def reference_selection(tree, overlap_dr):
    lepton_fields = ("PT", "Eta", "Phi", "Charge", "IsolationVar", "IsolationVarRhoCorr")
    electrons = records(tree, "Electron", *lepton_fields)
    muons = records(tree, "Muon", *lepton_fields)
    missing = records(tree, "MissingET", "MET", "Phi")
    large = records(tree, "LargeRJet", "PT", "Eta", "Phi", "Mass", "Constituents")
    small = records(tree, "SmallRJet", "PT", "Eta", "Phi", "Mass", "BTag")
    tracks = records(tree, "EFlowTrack", "fUniqueID", "PT", "Eta", "Phi", "Charge", "PID", "D0", "DZ", "ErrorD0", "ErrorDZ")
    photons = records(tree, "EFlowPhoton", "fUniqueID", "ET", "Eta", "Phi")
    neutrals = records(tree, "EFlowNeutralHadron", "fUniqueID", "ET", "Eta", "Phi")
    cutflow = [0, 0, 0, 0]
    accepted = {}
    for i, (els, mus, met, ljs, sjs) in enumerate(zip(electrons, muons, missing, large, small)):
        candidates = [(el, 11) for el in els if el["PT"] >= 30 and abs(el["Eta"]) <= 3]
        candidates += [(mu, 13) for mu in mus if mu["PT"] >= 30 and abs(mu["Eta"]) <= 3]
        choice = max(candidates, key=lambda pair: pair[0]["PT"], default=None)
        if choice is None:
            continue
        lepton, flavor = choice
        cutflow[0] += 1
        if met[0]["MET"] < 30:
            continue
        cutflow[1] += 1
        candidates = [jet for jet in ljs if jet["PT"] >= 250]
        if len(candidates) != 1 or abs(candidates[0]["Eta"]) > 3:
            continue
        cutflow[2] += 1
        large_jet = candidates[0]
        jets = sorted((jet for jet in sjs if dr(jet, large_jet) >= overlap_dr), key=lambda jet: -jet["PT"])
        if not any(jet["BTag"] & 1 and jet["PT"] > 30 and abs(jet["Eta"]) < 3 and 0.3 < dr(jet, lepton) < 3
                   for jet in jets):
            continue
        cutflow[3] += 1
        objects = {obj["fUniqueID"]: (obj, 1) for obj in tracks[i]}
        objects.update({obj["fUniqueID"]: (obj, 0) for obj in photons[i] + neutrals[i]})
        constituents = [(obj.get("PT", obj.get("ET")), obj, kind)
                        for obj, kind in (objects[uid] for uid in large_jet["Constituents"]["refs"])]
        constituents.sort(key=lambda row: -row[0])
        accepted[i] = (lepton, flavor, met[0], large_jet, constituents, jets)
    return cutflow, accepted


def check_array(output, row, name, values, capacity):
    assert output[name].shape[1] == capacity
    expected = list(values)[:capacity]
    np.testing.assert_array_equal(output[name][row], expected + [-999] * (capacity - len(expected)))


@pytest.mark.parametrize("nc,nj,overlap_dr", [(100, 10, 1.5), (3, 1, 1.5), (100, 10, 1.9)])
def test_debug_join_and_cutflow(debug_inputs, tmp_path, nc, nj, overlap_dr):
    mg, tag, reco, labels, large_r, small_r = debug_inputs
    # Full sample at defaults; one representative pair with capacities lowered
    # far enough that a b-tagged jet may pass the cut without being stored.
    runs = len(reco) if (nc, overlap_dr) == (100, 1.5) else 1
    env = analysis_env(ANALYSIS_OUTPUT_DIR=str(tmp_path), MAX_LARGE_R_CONSTITUENTS=str(nc),
                       MAX_SMALL_R_JETS=str(nj), JET_OVERLAP_DR=str(overlap_dr))
    subprocess.run([ANALYSIS / "run.sh", mg, tag, str(runs), "3"], cwd=REPO, env=env, check=True)

    for i in range(runs):
        stem = f"dataset_selected_{tag}_{i}"
        with uproot.open(f"{reco[i]}:Delphes") as tree:
            cutflow, expected = reference_selection(tree, overlap_dr)
        with uproot.open(f"{tmp_path / stem}.root:selected") as tree:
            output = tree.arrays(library="np")
        log = (tmp_path / "logs" / f"{stem}.log").read_text()
        for name, count in zip(("Passed lepton", "Passed MET", "Passed large-R jet",
                                "Passed leptonic b jet (accepted)"), cutflow):
            assert f"{name}: {count}\n" in log
        assert f"Accepted electrons: {sum(row[1] == 11 for row in expected.values())}\n" in log
        assert f"Accepted muons: {sum(row[1] == 13 for row in expected.values())}\n" in log
        assert "source_entry" not in output
        with uproot.open(f"{labels[i]}:labels") as tree:
            assert len(tree.keys()) == 50
            for name, values in tree.arrays(library="np").items():
                np.testing.assert_array_equal(output[name], values[list(expected)])
        for row, (lepton, flavor, met, large, constituents, jets) in enumerate(expected.values()):
            for name, field in [("pt", "PT"), ("eta", "Eta"), ("phi", "Phi"), ("charge", "Charge"),
                                ("isolation", "IsolationVar"), ("isolation_rho_corr", "IsolationVarRhoCorr")]:
                assert output[f"lepton_{name}"][row] == lepton[field]
            assert output["lepton_flavor"][row] == flavor
            assert output["met"][row] == met["MET"]
            assert output["large_r_jet_pt"][row] == large["PT"]
            for name, count, limit in [("constituents", len(constituents), nc), ("small_r_jets", len(jets), nj)]:
                assert output[f"num_{name}"][row] == count
                assert output[f"num_{name}_stored"][row] == min(count, limit)
                assert output[f"num_{name}_overflow"][row] == max(0, count - limit)
            check_array(output, row, "constituent_pt", (pt for pt, _, _ in constituents), nc)
            check_array(output, row, "constituent_type", (kind for _, _, kind in constituents), nc)
            for name, field in [("eta", "Eta"), ("phi", "Phi"), ("charge", "Charge"), ("pid", "PID"),
                                ("d0", "D0"), ("dz", "DZ"), ("d0_error", "ErrorD0"), ("dz_error", "ErrorDZ")]:
                check_array(output, row, f"constituent_{name}", (obj.get(field, -999) for _, obj, _ in constituents), nc)
            for name, field in [("pt", "PT"), ("eta", "Eta"), ("phi", "Phi"), ("mass", "Mass")]:
                check_array(output, row, f"small_r_jet_{name}", (jet[field] for jet in jets), nj)
            check_array(output, row, "small_r_jet_btag", (jet["BTag"] & 1 for jet in jets), nj)
        if nj == 1:
            assert (output["small_r_jet_btag"][:, 0] == 0).any()


def test_run_job_passes_unexported_config(debug_inputs, tmp_path):
    _, tag, _, _, _, _ = debug_inputs
    settings = {
        "ELECTRON_PT_MIN": "35", "ELECTRON_ABS_ETA_MAX": "2.4",
        "MUON_PT_MIN": "9999", "MUON_ABS_ETA_MAX": "2.8", "MET_MIN": "35",
        "LARGE_R_JET_PT_MIN": "260", "LARGE_R_JET_ABS_ETA_MAX": "2.9",
        "JET_OVERLAP_DR": "1.1", "LEPTONIC_B_JET_PT_MIN": "35",
        "LEPTONIC_B_JET_ABS_ETA_MAX": "2.8", "LEPTON_B_JET_DR_MIN": "0.4",
        "LEPTON_B_JET_DR_MAX": "2.9", "MAX_LARGE_R_CONSTITUENTS": "7",
        "MAX_SMALL_R_JETS": "4",
    }
    config = tmp_path / "job.config"
    config.write_text("source debug.config\nbypass_madgraph=true\nbypass_pythia=true\n"
                      "bypass_delphes=true\nbypass_analysis=false\nnum_runs=1\nmax_cpu_cores=1\n" +
                      "".join(f"{name}={value}\n" for name, value in settings.items()))
    subprocess.run([REPO / "run_job.sh", config], cwd=REPO,
                   env=analysis_env(ANALYSIS_OUTPUT_DIR=str(tmp_path)), check=True)
    stem = f"dataset_selected_{tag}_0"
    log = (tmp_path / "logs" / f"{stem}.log").read_text()
    actual = dict(line.split("=", 1) for line in log.splitlines()
                  if "=" in line and line.split("=", 1)[0] in settings)
    for name, value in settings.items():
        assert float(actual[name]) == float(value)
    with uproot.open(f"{tmp_path / stem}.root:selected") as tree:
        assert tree.num_entries > 0
        assert set(tree["lepton_flavor"].array().tolist()) == {11}
        assert tree["constituent_pt"].array(library="np").shape[1] == 7
        assert tree["small_r_jet_pt"].array(library="np").shape[1] == 4


@pytest.mark.parametrize("broken_input,message", [("count", "Entry-count mismatch"), ("schema", "incompatible MadGraph label")])
def test_invalid_labels_do_not_publish(debug_inputs, tmp_path, broken_input, message):
    _, _, reco, labels, large_r, small_r = debug_inputs
    with uproot.open(f"{labels[0]}:labels") as tree:
        data = tree.arrays(library="np")
    if broken_input == "count":
        data = {name: values[:-1] for name, values in data.items()}
    else:
        data.pop("W_E_topRest")
    bad = tmp_path / "bad.root"
    with uproot.recreate(bad) as file:
        file.mktree("labels", data)
    output = tmp_path / "selected.root"
    result = subprocess.run([ANALYSIS / "run_selections", reco[0], bad, output, str(large_r), str(small_r)],
                            env=analysis_env(), capture_output=True, text=True)
    assert result.returncode != 0
    assert message in result.stderr
    assert not output.exists()
    assert not list(tmp_path.glob("*.tmp.*"))


def make_boundary_sample(root, path, bad_sequence=False):
    # Exactly representable distances: overlap = 2, leptonic lower bound = 0.5.
    rows = [
        (1, 2, 31, False, False),      # tangency, retained
        (1, 1.875, 31, False, False),  # overlap, removed
        (1.5, 2, 31, False, False),    # b-jet lower DeltaR bound, rejected
        (-1, 2, 31, False, False),     # b-jet upper DeltaR bound, rejected
        (1, 2, 30, False, False),      # exclusive b-jet pT, rejected
        (1, 3, 31, False, False),      # exclusive b-jet eta, rejected
        (1, 2, 31, True, False),       # second hard jet outside eta acceptance, rejected
        (1, 0, 31, False, True),      # jet axes overlap across the phi seam
        (1, 2, 31, False, False),      # muon only
        (1, 2, 31, False, False),      # muon has higher pT
        (1, 2, 31, False, False),      # electron has higher pT
        (1, 2, 31, False, False),      # equal pT: electron wins
        (1, 2, 31, False, False),      # muon below threshold
        (0.5, 1.5, 31, False, False),  # exactly default overlap boundary
        (0.5, 1.375, 31, False, False), # just inside default overlap boundary
    ]
    classes = {"Event": "HepMCEvent", "Electron": "Electron", "Muon": "Muon", "MissingET": "MissingET",
               "LargeRJet": "Jet", "SmallRJet": "Jet", "EFlowTrack": "Track",
               "EFlowPhoton": "Tower", "EFlowNeutralHadron": "Tower"}
    with root.TFile(str(path), "RECREATE") as file:
        tree = root.TTree("Delphes", "Selection boundary cases")
        collections = {name: root.TClonesArray(cls) for name, cls in classes.items()}
        for name, arr in collections.items():
            assert tree.Branch(name, arr)

        def put(name, index=0, **fields):
            obj = collections[name].ConstructedAt(index)
            for key, value in fields.items():
                setattr(obj, key, value)

        for i, (electron_eta, b_eta, b_pt, extra_large, wraps) in enumerate(rows):
            for arr in collections.values():
                arr.Clear()
            put("Event", Number=i + int(bad_sequence and i == len(rows) - 1))
            if i not in (8, 12):
                put("Electron", PT=50 if i == 10 else 30, Eta=electron_eta, Phi=-3.1 if wraps else 0,
                    Charge=-1, IsolationVar=0.1, IsolationVarRhoCorr=0.1)
            if 8 <= i <= 12:
                muon_pt = {8: 30, 9: 50, 10: 30, 11: 30, 12: 29}[i]
                put("Muon", PT=muon_pt, Eta=electron_eta, Phi=0,
                    Charge=-1, IsolationVar=0.2, IsolationVarRhoCorr=0.2)
            put("MissingET", MET=30, Phi=0)
            put("LargeRJet", PT=250, Eta=0, Phi=3.1 if wraps else 0, Mass=100)
            if extra_large:
                put("LargeRJet", 1, PT=251, Eta=4, Phi=0, Mass=100)
            put("SmallRJet", PT=b_pt, Eta=b_eta, Phi=-3.1 if wraps else 0, Mass=10, BTag=1)
            tree.Fill()
        file.Write()
    return len(rows)


def test_boundaries_and_event_sequence(root_build, tmp_path):
    root = pytest.importorskip("ROOT")
    root.gInterpreter.AddIncludePath(str(DELPHES / "external"))
    assert root.gSystem.Load(str(DELPHES / "libDelphes.so")) >= 0
    reco = tmp_path / "delphes.root"
    entries = make_boundary_sample(root, reco)
    labels = tmp_path / "labels.root"
    names = [f"{obj}_{comp}_{frame}" for obj in ("top", "bottom", "down", "W")
             for comp in ("px", "py", "pz", "E") for frame in ("lab", "ttCM", "topRest")]
    names += ["costheta_down", "costheta_bottom"]
    with uproot.recreate(labels) as file:
        data = {name: np.zeros(entries, dtype=np.float32) for name in names}
        # Give each input a distinct test-only truth value to identify selected
        # entries without an event-index branch in the output.
        data["top_px_lab"] = np.arange(entries, dtype=np.float32)
        file.mktree("labels", data)

    output = tmp_path / "selected.root"
    cmd = [ANALYSIS / "run_selections", reco, labels, output, "1.5", "0.5"]
    env = analysis_env(LEPTON_B_JET_DR_MIN="0.5", JET_OVERLAP_DR="2.0")
    subprocess.run(cmd, env=env, check=True)
    with uproot.open(f"{output}:selected") as tree:
        assert "source_entry" not in tree
        documented = []
        rows = re.findall(
            r"^\| `([A-Za-z0-9_,{}]+)` \| `(Float_t|Int_t)(?:\[(Nc|Nj)\])?` \|",
            (ANALYSIS / "README.md").read_text(), flags=re.MULTILINE,
        )
        for pattern, dtype, capacity in rows:
            grouped = re.fullmatch(r"(\w+)_\{([^{}]+)\}_\{([^{}]+)\}", pattern)
            if grouped:
                particle, components, frames = grouped.groups()
                documented.extend((f"{particle}_{component}_{frame}", dtype, capacity)
                                  for component in components.split(",") for frame in frames.split(","))
            else:
                assert "{" not in pattern and "}" not in pattern
                documented.append((pattern, dtype, capacity))
        assert len(documented) == 85
        assert set(name for name, _, _ in documented) == set(tree.keys())
        for name, dtype, capacity in documented:
            values = tree[name].array(library="np")
            assert values.dtype == {"Float_t": np.dtype("float32"), "Int_t": np.dtype("int32")}[dtype]
            shape = {"": (), "Nc": (100,), "Nj": (10,)}[capacity]
            assert values.shape == (tree.num_entries, *shape)
        assert tree["top_px_lab"].array().tolist() == [0, 8, 9, 10, 11]
        assert tree["lepton_flavor"].array().tolist() == [11, 13, 13, 11, 11]
        assert tree["lepton_pt"].array().tolist() == [30, 30, 50, 50, 30]
        assert tree["num_constituents"].array().tolist() == [0] * 5
        assert tree["constituent_pt"].array().tolist() == [[-999] * 100] * 5

    default_output = tmp_path / "selected_default.root"
    subprocess.run([ANALYSIS / "run_selections", reco, labels, default_output, "1.5", "0.5"],
                   env=analysis_env(LEPTON_B_JET_DR_MIN="0.5"), check=True)
    with uproot.open(f"{default_output}:selected") as tree:
        assert tree["top_px_lab"].array().tolist() == [0, 1, 8, 9, 10, 11, 13]

    # Flavor thresholds can be changed independently. At a pT tie the first
    # eligible electron wins; raising its threshold leaves the muon eligible.
    threshold_output = tmp_path / "selected_threshold.root"
    subprocess.run([ANALYSIS / "run_selections", reco, labels, threshold_output, "1.5", "0.5"],
                   env=analysis_env(LEPTON_B_JET_DR_MIN="0.5", JET_OVERLAP_DR="2.0",
                                    ELECTRON_PT_MIN="40"), check=True)
    with uproot.open(f"{threshold_output}:selected") as tree:
        assert tree["top_px_lab"].array().tolist() == [8, 9, 10, 11]
        assert tree["lepton_flavor"].array().tolist() == [13, 13, 11, 13]

    subprocess.run([ANALYSIS / "run_selections", reco, labels, threshold_output, "1.5", "0.5"],
                   env=analysis_env(LEPTON_B_JET_DR_MIN="0.5", JET_OVERLAP_DR="2.0",
                                    MUON_PT_MIN="40"), check=True)
    with uproot.open(f"{threshold_output}:selected") as tree:
        assert tree["top_px_lab"].array().tolist() == [0, 9, 10, 11]
        assert tree["lepton_flavor"].array().tolist() == [11, 13, 11, 11]

    previous = output.read_bytes()
    make_boundary_sample(root, reco, bad_sequence=True)
    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    assert result.returncode != 0
    assert "Non-sequential Delphes Event.Number" in result.stderr
    assert output.read_bytes() == previous
    assert not list(tmp_path.glob("*.tmp.*"))
