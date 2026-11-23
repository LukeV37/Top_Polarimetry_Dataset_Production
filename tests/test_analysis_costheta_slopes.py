"""Fit and validate truth cos(theta) in a combined selected ROOT tree.

Run with --file and --polarization; optional input files check the merge and
compare against the original, preselection MadGraph labels. See analysis/README.md.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import uproot

from test_madgraph_costheta_slopes import EXPECTED_SLOPES, _fit_histogram_beta


ANGLES = ("costheta_down", "costheta_bottom")
VECTOR_BRANCHES = tuple(
    f"{particle}_{component}_{frame}"
    for particle, frame in (("top", "ttCM"), ("down", "topRest"), ("bottom", "topRest"))
    for component in ("px", "py", "pz")
)


def read_angles(path: Path, tree_name: str, check_vectors: bool = False,
                allow_empty: bool = False) -> dict[str, np.ndarray]:
    branches = ANGLES + VECTOR_BRANCHES if check_vectors else ANGLES
    with uproot.open(path) as root_file:
        arrays = root_file[tree_name].arrays(branches, library="np")
        values = {name: np.asarray(arrays[name], dtype=np.float64) for name in branches}

    if not allow_empty and not len(values[ANGLES[0]]):
        raise ValueError(f"No events in {path}:{tree_name}")
    for name, array in values.items():
        if array.ndim != 1 or not np.isfinite(array).all():
            raise ValueError(f"Invalid or non-finite {name} in {path}")
    for name in ANGLES:
        if np.any(np.abs(values[name]) > 1.0):
            raise ValueError(f"{name} outside [-1, 1] in {path}")

    if check_vectors:
        top = np.column_stack([values[f"top_{c}_ttCM"] for c in ("px", "py", "pz")])
        for particle in ("down", "bottom"):
            daughter = np.column_stack([values[f"{particle}_{c}_topRest"] for c in ("px", "py", "pz")])
            norms = np.linalg.norm(top, axis=1) * np.linalg.norm(daughter, axis=1)
            if np.any(norms == 0) or not np.isfinite(norms).all():
                raise ValueError(f"Undefined {particle} or top direction in {path}")
            reconstructed = np.sum(top * daughter, axis=1) / norms
            error = np.max(np.abs(reconstructed - values[f"costheta_{particle}"]))
            if error > 1e-5:
                raise ValueError(f"costheta_{particle} disagrees with stored vectors in {path}: max error {error:g}")
    return {name: values[name] for name in ANGLES}


def read_many(paths: list[Path], tree_name: str) -> dict[str, np.ndarray]:
    if not paths:
        raise ValueError("No input files supplied")
    chunks = [read_angles(path, tree_name, allow_empty=True) for path in paths]
    values = {name: np.concatenate([chunk[name] for chunk in chunks]) for name in ANGLES}
    if not len(values[ANGLES[0]]):
        raise ValueError(f"No events across the input {tree_name} files")
    return values


def angle_pairs(values: dict[str, np.ndarray]) -> np.ndarray:
    # Keep the two labels from each event together; do not assume hadd order.
    return np.rec.fromarrays([values[name] for name in ANGLES], names=ANGLES)


def check_merge(selected: dict[str, np.ndarray], parts: dict[str, np.ndarray]) -> None:
    if not np.array_equal(np.sort(angle_pairs(selected)), np.sort(angle_pairs(parts))):
        raise ValueError("Combined angle pairs do not match the per-run selected files")


def check_label_subset(selected: dict[str, np.ndarray], labels: dict[str, np.ndarray]) -> None:
    selected_pairs, selected_counts = np.unique(angle_pairs(selected), return_counts=True)
    label_pairs, label_counts = np.unique(angle_pairs(labels), return_counts=True)
    indices = np.searchsorted(label_pairs, selected_pairs)
    if np.any(indices == len(label_pairs)):
        raise ValueError("Selected angle pairs are not a subset of the MadGraph labels")
    if np.any(label_pairs[indices] != selected_pairs) or np.any(selected_counts > label_counts[indices]):
        raise ValueError("Selected angle pairs are not a subset of the MadGraph labels")


def print_fits(values: dict[str, np.ndarray], polarization: str, bins: int, tolerance: float) -> bool:
    all_match = True
    for name in ANGLES:
        beta, intercept = _fit_histogram_beta(values[name], bins=bins)
        reference = EXPECTED_SLOPES[polarization][name]
        matches = abs(beta - reference) <= tolerance
        all_match &= matches
        print(f"  {name}: beta={beta:+.4f}, intercept={intercept:.4f}; "
              f"preselection reference={reference:+.2f} ± {tolerance:.2f} "
              f"({'within' if matches else 'outside'} tolerance)")
    return all_match


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", required=True, type=Path, help="Combined ROOT file with the selected tree.")
    parser.add_argument("--polarization", required=True, choices=("U", "L", "R"))
    parser.add_argument("--run-files", nargs="+", type=Path, help="Per-run selected ROOT files to compare with the merge.")
    parser.add_argument("--label-files", nargs="+", type=Path, help="Original MadGraph label ROOT files to compare.")
    parser.add_argument("--bins", type=int, default=60, help="Histogram bins (default: 60, as in the MadGraph test).")
    parser.add_argument("--reference-tolerance", type=float, default=0.10,
                        help="Comparison to preselection reference beta (default: 0.10).")
    parser.add_argument("--strict-reference", action="store_true",
                        help="Fail if selected beta differs from the preselection reference by more than the tolerance.")
    args = parser.parse_args()
    if args.bins < 2 or not np.isfinite(args.reference_tolerance) or args.reference_tolerance < 0:
        parser.error("--bins must be at least 2 and --reference-tolerance must be finite and nonnegative")

    selected = read_angles(args.file, "selected", check_vectors=True)
    print(f"{args.polarization}: {args.file}\nSelected entries: {len(selected[ANGLES[0]])}")
    print("Stored cos(theta) agrees with the stored frame vectors (absolute tolerance 1e-5).")
    print("Selected fits: f(x) = intercept + 0.5 * beta * x")
    matches = print_fits(selected, args.polarization, args.bins, args.reference_tolerance)

    if args.run_files:
        parts = read_many(args.run_files, "selected")
        check_merge(selected, parts)
        print(f"Merge check: angle pairs match {len(args.run_files)} per-run files "
              f"({len(parts[ANGLES[0]])} entries), independent of order.")
    if args.label_files:
        labels = read_many(args.label_files, "labels")
        check_label_subset(selected, labels)
        print(f"MadGraph labels: {len(labels[ANGLES[0]])} entries from {len(args.label_files)} files; "
              "all selected angle pairs are present.")
        print("Preselection fits:")
        print_fits(labels, args.polarization, args.bins, args.reference_tolerance)

    print("Reference values describe preselection truth; detector/analysis selection can change the shape.")
    if args.strict_reference and not matches:
        print("Selected beta is outside the requested preselection-reference tolerance.", file=sys.stderr)
        return 1
    return 0


def test_angle_pair_checks_preserve_event_pairs_and_multiplicity() -> None:
    selected = {ANGLES[0]: np.array([0.1, 0.2, 0.1]), ANGLES[1]: np.array([0.3, 0.4, 0.3])}
    reordered = {name: values[[1, 0, 2]] for name, values in selected.items()}
    check_merge(selected, reordered)
    check_label_subset(selected, reordered)
    wrong_pairs = {ANGLES[0]: selected[ANGLES[0]], ANGLES[1]: np.array([0.4, 0.3, 0.3])}
    with pytest.raises(ValueError, match="per-run selected"):
        check_merge(selected, wrong_pairs)
    with pytest.raises(ValueError, match="subset"):
        check_label_subset(selected, wrong_pairs)
    insufficient_copies = {name: values[[0, 1, 1]] for name, values in selected.items()}
    with pytest.raises(ValueError, match="subset"):
        check_label_subset(selected, insufficient_copies)


def test_read_angles_validates_frames_and_values(tmp_path: Path) -> None:
    path = tmp_path / "selected.root"
    values = {name: np.array([0.0, 0.0], dtype=np.float32) for name in ANGLES + VECTOR_BRANCHES}
    values["top_px_ttCM"][:] = 1
    values["down_py_topRest"][:] = 1
    values["bottom_px_topRest"][:] = 1
    values["costheta_bottom"][:] = 1
    with uproot.recreate(path) as root_file:
        root_file.mktree("selected", values)
    assert len(read_angles(path, "selected", check_vectors=True)[ANGLES[0]]) == 2
    values["costheta_down"][0] = 0.5
    with uproot.recreate(path) as root_file:
        root_file.mktree("selected", values)
    with pytest.raises(ValueError, match="disagrees with stored vectors"):
        read_angles(path, "selected", check_vectors=True)
    values["costheta_down"][0] = np.nan
    with uproot.recreate(path) as root_file:
        root_file.mktree("selected", values)
    with pytest.raises(ValueError, match="non-finite"):
        read_angles(path, "selected", check_vectors=True)
    values["costheta_down"][0] = 1.01
    with uproot.recreate(path) as root_file:
        root_file.mktree("selected", values)
    with pytest.raises(ValueError, match=r"outside \[-1, 1\]"):
        read_angles(path, "selected", check_vectors=True)
    values["costheta_down"][0] = 0
    values["top_px_ttCM"][0] = 0
    with uproot.recreate(path) as root_file:
        root_file.mktree("selected", values)
    with pytest.raises(ValueError, match="Undefined"):
        read_angles(path, "selected", check_vectors=True)


def test_cli_fits_and_checks_inputs(tmp_path: Path) -> None:
    # A flat, exactly populated histogram, with vectors consistent with angles.
    angles = np.linspace(-59 / 60, 59 / 60, 60, dtype=np.float32)
    values = {name: np.zeros(60, dtype=np.float32) for name in ANGLES + VECTOR_BRANCHES}
    values["top_px_ttCM"][:] = 1
    for particle in ("down", "bottom"):
        values[f"costheta_{particle}"] = angles.copy()
        values[f"{particle}_px_topRest"] = angles.copy()
        values[f"{particle}_py_topRest"] = np.sqrt(1 - angles**2)

    combined, part, labels = (tmp_path / name for name in ("combined.root", "part.root", "labels.root"))
    for path, tree, data in (
        (combined, "selected", values),
        (part, "selected", {name: values[name][::-1] for name in ANGLES}),
        (labels, "labels", {name: values[name] for name in ANGLES}),
    ):
        with uproot.recreate(path) as root_file:
            root_file.mktree(tree, data)
    command = [sys.executable, str(Path(__file__).resolve()), "--file", str(combined),
               "--run-files", str(part), "--label-files", str(labels)]
    result = subprocess.run(command + ["--polarization", "U", "--strict-reference"],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "Selected entries: 60" in result.stdout
    assert "Merge check: angle pairs match 1 per-run files" in result.stdout
    assert "all selected angle pairs are present" in result.stdout
    result = subprocess.run(command + ["--polarization", "R", "--strict-reference"],
                            capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "outside the requested preselection-reference tolerance" in result.stderr


if __name__ == "__main__":
    raise SystemExit(main())
