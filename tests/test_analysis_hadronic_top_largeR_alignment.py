"""Validate the selected large-R jet against the joined MadGraph top direction.

Run directly with --analysis-dir, or use pytest with the workspace in
debug.config. TP_ANALYSIS_DIR, TP_CONFIG, TP_LARGE_R_RADIUS, and
TP_MIN_MATCH_FRACTION customize the pytest workflow.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
import uproot


REPO = Path(__file__).resolve().parents[1]
FILE_PATTERN = re.compile(r"dataset_selected_(.+)_(\d+)\.root\Z")
RADIUS_PATTERN = re.compile(r"^R_LARGE=(\S+)\s*$", re.MULTILINE)
BRANCHES = ("top_px_lab", "top_py_lab", "top_pz_lab", "large_r_jet_eta", "large_r_jet_phi")


@dataclass(frozen=True)
class AlignmentSummary:
    analysis_dir: Path
    jet_radius: float
    file_count: int
    total_events: int
    matched_events: int
    match_fraction: float
    p95_delta_r: float
    max_delta_r: float


def _selected_files(directory: Path) -> list[Path]:
    indexed = {}
    tags = set()
    for path in directory.glob("dataset_selected_*.root"):
        match = FILE_PATTERN.fullmatch(path.name)
        if not match:
            raise ValueError(f"Cannot parse selected file name: {path}")
        tag, index = match.group(1), int(match.group(2))
        if index in indexed:
            raise ValueError(f"Duplicate run index {index}: {indexed[index]} and {path}")
        indexed[index] = path
        tags.add(tag)
    if not indexed:
        raise FileNotFoundError(f"No selected ROOT files found in {directory}")
    if len(tags) != 1 or sorted(indexed) != list(range(len(indexed))):
        raise ValueError(f"Expected one tag and contiguous run indices from 0 in {directory}")
    return [indexed[index] for index in sorted(indexed)]


def _positive_radius(value: float) -> float:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"Jet radius must be finite and positive, got {value}")
    return value


def _radius_from_logs(files: list[Path]) -> float:
    radii = []
    for path in files:
        log = path.parent / "logs" / f"{path.stem}.log"
        if not log.is_file():
            raise FileNotFoundError(f"Missing analysis log {log}; pass --jet-radius if logs are unavailable")
        matches = RADIUS_PATTERN.findall(log.read_text())
        if len(matches) != 1:
            raise ValueError(f"Expected one R_LARGE in {log}; pass --jet-radius if logs are unavailable")
        radii.append(_positive_radius(float(matches[0])))
    if len(set(radii)) != 1:
        raise ValueError(f"Inconsistent R_LARGE values in analysis logs: {radii}")
    return radii[0]


def summarize_alignment(analysis_dir: Path, jet_radius: float | None = None) -> AlignmentSummary:
    analysis_dir = analysis_dir.resolve()
    files = _selected_files(analysis_dir)
    radius = _radius_from_logs(files) if jet_radius is None else _positive_radius(jet_radius)
    distances = []
    for path in files:
        with uproot.open(path) as root_file:
            tree = root_file["selected"]
            values = tree.arrays(BRANCHES, library="np")
            px, py, pz, jet_eta, jet_phi = (np.asarray(values[name], dtype=np.float64) for name in BRANCHES)
            if any(not np.all(np.isfinite(array)) for array in (px, py, pz, jet_eta, jet_phi)):
                raise ValueError(f"Non-finite top or jet direction in {path}")
            pt = np.hypot(px, py)
            if np.any(pt == 0):
                raise ValueError(f"Hadronic top has undefined pseudorapidity in {path}")
            top_eta = np.arcsinh(pz / pt)
            top_phi = np.arctan2(py, px)
            delta_phi = np.arctan2(np.sin(top_phi - jet_phi), np.cos(top_phi - jet_phi))
            delta_r = np.hypot(top_eta - jet_eta, delta_phi)
            if not np.all(np.isfinite(delta_r)):
                raise ValueError(f"Non-finite top/jet ΔR in {path}")
            distances.append(delta_r)

    delta_r = np.concatenate(distances)
    total = int(delta_r.size)
    if not total:
        raise ValueError(f"No selected events found in {analysis_dir}")
    matched = int(np.count_nonzero(delta_r <= radius))
    return AlignmentSummary(analysis_dir, radius, len(files), total, matched, matched / total,
                            float(np.quantile(delta_r, 0.95)), float(np.max(delta_r)))


def assert_alignment_tolerance(summary: AlignmentSummary, min_match_fraction: float) -> None:
    if not 0 <= min_match_fraction <= 1:
        raise ValueError(f"min_match_fraction must be between 0 and 1, got {min_match_fraction}")
    if summary.match_fraction < min_match_fraction:
        raise AssertionError(
            f"Hadronic-top / selected-large-R-jet alignment fell below tolerance: "
            f"{summary.matched_events} / {summary.total_events} matched within ΔR <= {summary.jet_radius:g} "
            f"({summary.match_fraction:.4%}), below {min_match_fraction:.4%}."
        )


def _format_summary(summary: AlignmentSummary) -> str:
    return "\n".join([
        f"analysis dir: {summary.analysis_dir}",
        f"files checked: {summary.file_count}",
        f"Large-R radius used for match test: {summary.jet_radius:.4f}",
        f"selected events: {summary.total_events}",
        f"matched events with ΔR(top, selected Large-R jet) <= R: {summary.matched_events}",
        f"unmatched selected events: {summary.total_events - summary.matched_events}",
        f"match fraction: {summary.match_fraction:.4%}",
        f"95th percentile ΔR: {summary.p95_delta_r:.4f}",
        f"max ΔR: {summary.max_delta_r:.4f}",
    ])


def _debug_analysis_dir() -> Path:
    config = Path(os.getenv("TP_CONFIG", REPO / "debug.config")).resolve()
    tag = subprocess.check_output(
        ["bash", "-c", 'set -e; source "$1"; printf "%s" "$DELPHES_tag"', "analysis-test-config", str(config)],
        cwd=REPO, text=True,
    )
    if not tag:
        raise ValueError(f"Expected DELPHES_tag in {config}")
    return REPO / "analysis" / f"WS_{tag}"


def test_selected_hadronic_top_large_r_alignment() -> None:
    override = os.getenv("TP_ANALYSIS_DIR")
    analysis_dir = Path(override) if override else _debug_analysis_dir()
    if not override and not list(analysis_dir.glob("dataset_selected_*.root")):
        pytest.skip("Requires existing selected outputs from debug.config")
    radius = os.getenv("TP_LARGE_R_RADIUS")
    summary = summarize_alignment(analysis_dir, float(radius) if radius is not None else None)
    print(_format_summary(summary))
    assert_alignment_tolerance(summary, float(os.getenv("TP_MIN_MATCH_FRACTION", "0.90")))


def _write_sample(directory: Path, index: int, rows: list[tuple[float, ...]], radius: float = 1.0) -> None:
    stem = f"dataset_selected_sample_{index}"
    data = np.asarray(rows, dtype=np.float64).reshape(-1, len(BRANCHES))
    with uproot.recreate(directory / f"{stem}.root") as root_file:
        root_file.mktree("selected", {name: data[:, i] for i, name in enumerate(BRANCHES)})
    (directory / "logs").mkdir(exist_ok=True)
    (directory / "logs" / f"{stem}.log").write_text(f"R_LARGE={radius}\n")


def test_alignment_wrap_boundary_and_tolerance(tmp_path: Path) -> None:
    # First event crosses the phi seam; the second is exactly at ΔR = R.
    _write_sample(tmp_path, 0, [(-1, 0.01, 0, 0, -math.pi + 0.01), (1, 0, 0, 1, 0)])
    _write_sample(tmp_path, 1, [(1, 0, 0, 1.01, 0), (1, 0, 0, 0, math.pi / 2)])
    summary = summarize_alignment(tmp_path)
    assert (summary.file_count, summary.total_events, summary.matched_events) == (2, 4, 2)
    assert summary.match_fraction == 0.5
    assert summary.max_delta_r == pytest.approx(math.pi / 2)
    assert_alignment_tolerance(summary, 0.5)
    with pytest.raises(AssertionError, match="2 / 4 matched"):
        assert_alignment_tolerance(summary, 0.5001)
    with pytest.raises(ValueError, match="between 0 and 1"):
        assert_alignment_tolerance(summary, float("nan"))

    command = [sys.executable, str(Path(__file__).resolve()), "--analysis-dir", str(tmp_path)]
    result = subprocess.run(command + ["--min-match-fraction", "0.5"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "match fraction: 50.0000%" in result.stdout
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert "below tolerance" in result.stderr


@pytest.mark.parametrize("rows,message", [
    ([], "No selected events"),
    ([(0, 0, 1, 0, 0)], "undefined pseudorapidity"),
    ([(1, 0, 0, float("nan"), 0)], "Non-finite"),
])
def test_alignment_rejects_unusable_events(tmp_path: Path, rows, message: str) -> None:
    _write_sample(tmp_path, 0, rows)
    with pytest.raises(ValueError, match=message):
        summarize_alignment(tmp_path)


def test_alignment_radius_and_file_validation(tmp_path: Path) -> None:
    _write_sample(tmp_path, 0, [(1, 0, 0, 0, 0)])
    _write_sample(tmp_path, 1, [(1, 0, 0, 0, 0)], radius=1.5)
    with pytest.raises(ValueError, match="Inconsistent R_LARGE"):
        summarize_alignment(tmp_path)
    assert summarize_alignment(tmp_path, 1.0).matched_events == 2
    with pytest.raises(ValueError, match="finite and positive"):
        summarize_alignment(tmp_path, float("nan"))
    (tmp_path / "logs" / "dataset_selected_sample_1.log").unlink()
    with pytest.raises(FileNotFoundError, match="pass --jet-radius"):
        summarize_alignment(tmp_path)
    (tmp_path / "dataset_selected_sample_1.root").rename(tmp_path / "dataset_selected_sample_2.root")
    with pytest.raises(ValueError, match="contiguous"):
        summarize_alignment(tmp_path, 1.0)


def _cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate selected jet / MadGraph hadronic-top alignment.")
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--jet-radius", type=float, default=None,
                        help="Match radius; by default read R_LARGE from each analysis log.")
    parser.add_argument("--min-match-fraction", type=float, default=0.90,
                        help="Minimum fraction of selected events within ΔR <= radius (default: 0.90).")
    return parser


def main() -> int:
    args = _cli_parser().parse_args()
    summary = summarize_alignment(args.analysis_dir, args.jet_radius)
    print(_format_summary(summary))
    assert_alignment_tolerance(summary, args.min_match_fraction)
    print("Alignment test passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
