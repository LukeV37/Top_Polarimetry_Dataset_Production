"""Validate the selected large-R jet against the joined MadGraph top direction.

Run directly with --analysis-dir or --file, or use pytest with the workspace
in debug.config. TP_ANALYSIS_DIR, TP_CONFIG, TP_LARGE_R_RADIUS, and
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
DELTA_R_BINS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, np.inf)


@dataclass(frozen=True)
class AlignmentSummary:
    source: Path
    jet_radius: float
    file_count: int
    total_events: int
    matched_events: int
    match_fraction: float
    min_delta_r: float
    mean_delta_r: float
    std_delta_r: float
    median_delta_r: float
    p90_delta_r: float
    p95_delta_r: float
    p99_delta_r: float
    max_delta_r: float
    histogram: tuple[int, ...]


def _selected_files(directory: Path) -> list[Path]:
    indexed = {}
    tags = set()
    # An optional hadd output in the usual WS_<tag> directory must not be
    # counted in addition to the individual runs.
    combined = f"dataset_selected_{directory.name[3:]}.root" if directory.name.startswith("WS_") else None
    for path in directory.glob("dataset_selected_*.root"):
        if path.name == combined:
            continue
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
    return _summarize_files(analysis_dir, files, radius)


def summarize_file(path: Path, jet_radius: float) -> AlignmentSummary:
    path = path.resolve()
    return _summarize_files(path, [path], _positive_radius(jet_radius))


def _summarize_files(source: Path, files: list[Path], radius: float) -> AlignmentSummary:
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
        raise ValueError(f"No selected events found in {source}")
    matched = int(np.count_nonzero(delta_r <= radius))
    counts, _ = np.histogram(delta_r, bins=DELTA_R_BINS)
    return AlignmentSummary(
        source=source, jet_radius=radius, file_count=len(files), total_events=total,
        matched_events=matched, match_fraction=matched / total,
        min_delta_r=float(np.min(delta_r)), mean_delta_r=float(np.mean(delta_r)),
        std_delta_r=float(np.std(delta_r)), median_delta_r=float(np.median(delta_r)),
        p90_delta_r=float(np.quantile(delta_r, 0.90)),
        p95_delta_r=float(np.quantile(delta_r, 0.95)),
        p99_delta_r=float(np.quantile(delta_r, 0.99)), max_delta_r=float(np.max(delta_r)),
        histogram=tuple(int(count) for count in counts),
    )


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
    lines = [
        f"source: {summary.source}",
        f"files checked: {summary.file_count}",
        f"Large-R radius used for match test: {summary.jet_radius:.4f}",
        f"selected events: {summary.total_events}",
        f"matched events with ΔR(top, selected Large-R jet) <= R: {summary.matched_events}",
        f"unmatched selected events: {summary.total_events - summary.matched_events}",
        f"match fraction: {summary.match_fraction:.4%}",
        f"min ΔR: {summary.min_delta_r:.4f}",
        f"mean ΔR: {summary.mean_delta_r:.4f}",
        f"population std ΔR: {summary.std_delta_r:.4f}",
        f"median ΔR: {summary.median_delta_r:.4f}",
        f"90th percentile ΔR: {summary.p90_delta_r:.4f}",
        f"95th percentile ΔR: {summary.p95_delta_r:.4f}",
        f"99th percentile ΔR: {summary.p99_delta_r:.4f}",
        f"max ΔR: {summary.max_delta_r:.4f}",
        "ΔR histogram (lower edge inclusive, upper edge exclusive):",
    ]
    for lower, upper, count in zip(DELTA_R_BINS, DELTA_R_BINS[1:], summary.histogram):
        lines.append(f"  [{lower:g}, {upper:g}): {count} ({count / summary.total_events:.2%})")
    return "\n".join(lines)


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
    distances = np.array([math.atan(0.01) + 0.01, 1.0, 1.01, math.pi / 2])
    assert summary.min_delta_r == pytest.approx(distances.min())
    assert summary.mean_delta_r == pytest.approx(distances.mean())
    assert summary.std_delta_r == pytest.approx(distances.std())
    assert summary.median_delta_r == pytest.approx(np.median(distances))
    assert summary.p90_delta_r == pytest.approx(np.quantile(distances, 0.9))
    assert summary.p95_delta_r == pytest.approx(np.quantile(distances, 0.95))
    assert summary.p99_delta_r == pytest.approx(np.quantile(distances, 0.99))
    assert summary.max_delta_r == pytest.approx(math.pi / 2)
    assert summary.histogram == (1, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0)
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


def test_combined_file_and_directory_do_not_double_count(tmp_path: Path) -> None:
    directory = tmp_path / "WS_sample"
    directory.mkdir()
    _write_sample(directory, 0, [(1, 0, 0, 0, 0)])
    _write_sample(directory, 1, [(1, 0, 0, 2, 0)])
    combined = directory / "dataset_selected_sample.root"
    rows = np.array([(1, 0, 0, 0, 0), (1, 0, 0, 2, 0)], dtype=np.float64)
    with uproot.recreate(combined) as root_file:
        root_file.mktree("selected", {name: rows[:, i] for i, name in enumerate(BRANCHES)})

    from_directory = summarize_alignment(directory)
    from_file = summarize_file(combined, 1.0)
    assert from_directory.file_count == 2
    assert from_file.file_count == 1
    assert from_directory.total_events == from_file.total_events == 2
    assert from_directory.histogram == from_file.histogram
    assert from_directory.match_fraction == from_file.match_fraction == 0.5
    command = [sys.executable, str(Path(__file__).resolve()), "--file", str(combined)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert "--jet-radius is required" in result.stderr
    result = subprocess.run(command + ["--jet-radius", "1", "--min-match-fraction", "0.5"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "selected events: 2" in result.stdout
    assert "[2, 3): 1 (50.00%)" in result.stdout


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
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--analysis-dir", type=Path, help="Directory containing per-run selected files.")
    group.add_argument("--file", type=Path, help="Single selected ROOT file, such as a combined output.")
    parser.add_argument("--jet-radius", type=float, default=None,
                        help="Match radius; required with --file, otherwise read R_LARGE from analysis logs.")
    parser.add_argument("--min-match-fraction", type=float, default=0.90,
                        help="Minimum fraction of selected events within ΔR <= radius (default: 0.90).")
    return parser


def main() -> int:
    parser = _cli_parser()
    args = parser.parse_args()
    if args.file is not None:
        if args.jet_radius is None:
            parser.error("--jet-radius is required with --file")
        summary = summarize_file(args.file, args.jet_radius)
    else:
        summary = summarize_alignment(args.analysis_dir, args.jet_radius)
    print(_format_summary(summary))
    assert_alignment_tolerance(summary, args.min_match_fraction)
    print("Alignment test passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
