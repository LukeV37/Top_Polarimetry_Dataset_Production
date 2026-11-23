"""Check that the hadronic-top truth direction is captured by the leading Large-R jet.

This file supports two workflows:

1. Direct CLI validation over a pair of output directories:

   python3 tests/test_delphes_hadronic_top_largeR_alignment.py \
       --labels-dir madgraph/pp_tt_semi_full_<MG_TAG> \
       --delphes-dir delphes/WS_<DELPHES_TAG>

2. Local pytest execution, using the output tags in ``debug.config`` by default:

   python3 -m pytest -q -s tests/test_delphes_hadronic_top_largeR_alignment.py

   Set ``TP_LABELS_DIR`` and ``TP_DELPHES_DIR`` to override those defaults.

The main assertion is conditional on a leading Large-R jet existing in the
event. Events without any Large-R jet are reported separately and can also be
bounded with an optional tolerance.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import awkward as ak
import numpy as np
import pytest
import uproot


LARGE_R_RADIUS_PATTERN = re.compile(r"^\s*set\s+LargeRJetRadius\s+([0-9eE+\-.]+)\s*$")


@dataclass(frozen=True)
class AlignmentSummary:
    labels_dir: Path
    delphes_dir: Path
    jet_radius: float
    file_count: int
    total_events: int
    events_with_large_r_jet: int
    events_without_large_r_jet: int
    matched_events: int
    unmatched_events_with_large_r_jet: int
    inclusive_match_fraction: float
    conditional_match_fraction: float
    conditional_miss_fraction: float
    max_delta_r_with_large_r_jet: float
    p95_delta_r_with_large_r_jet: float


def _run_index(path: Path, prefix: str) -> int:
    stem = path.stem
    if not stem.startswith(prefix):
        raise ValueError(f"Unexpected file name '{path.name}'; expected prefix '{prefix}'")
    run_token = stem.rsplit("_", 1)[-1]
    if not run_token.isdigit():
        raise ValueError(f"Could not parse run index from '{path.name}'")
    return int(run_token)


def _indexed_root_files(directory: Path, prefix: str) -> dict[int, Path]:
    files = sorted(directory.glob(f"{prefix}*.root"))
    indexed = {_run_index(path, prefix): path for path in files}
    if not indexed:
        raise FileNotFoundError(f"No ROOT files matching '{prefix}*.root' found in {directory}")
    expected = list(range(len(indexed)))
    found = sorted(indexed)
    if found != expected:
        raise ValueError(
            f"Expected contiguous run indices in {directory}; found {found}, expected {expected}"
        )
    return indexed


def _infer_large_r_radius(delphes_dir: Path) -> float:
    config_dir = delphes_dir / "config"
    card_paths = sorted(config_dir.glob("delphes_card_*.tcl"))
    if len(card_paths) != 1:
        raise ValueError(
            "Could not infer the Large-R radius automatically. "
            f"Expected exactly one rendered Delphes card in {config_dir}, found {len(card_paths)}. "
            "Pass --jet-radius explicitly."
        )

    for line in card_paths[0].read_text().splitlines():
        match = LARGE_R_RADIUS_PATTERN.match(line)
        if match:
            return float(match.group(1))

    raise ValueError(f"LargeRJetRadius not found in {card_paths[0]}")


def _eta_from_momentum(px: np.ndarray, py: np.ndarray, pz: np.ndarray) -> np.ndarray:
    pt = np.hypot(px, py)
    eta = np.full_like(pt, np.nan, dtype=np.float64)
    mask = pt > 0.0
    eta[mask] = np.arcsinh(pz[mask] / pt[mask])
    return eta


def _delta_phi(phi1: np.ndarray, phi2: np.ndarray) -> np.ndarray:
    return np.arctan2(np.sin(phi1 - phi2), np.cos(phi1 - phi2))


def _iter_alignment_delta_r(
    labels_files: dict[int, Path], delphes_files: dict[int, Path]
) -> Iterable[np.ndarray]:
    label_indices = sorted(labels_files)
    delphes_indices = sorted(delphes_files)
    if label_indices != delphes_indices:
        raise ValueError(
            "Run-index mismatch between labels and Delphes outputs: "
            f"labels={label_indices}, delphes={delphes_indices}"
        )

    for run_index in label_indices:
        labels_root = labels_files[run_index]
        delphes_root = delphes_files[run_index]

        with uproot.open(labels_root) as label_file, uproot.open(delphes_root) as delphes_file:
            labels_tree = label_file["labels"]
            delphes_tree = delphes_file["Delphes"]

            label_vectors = labels_tree.arrays(
                ["top_px_lab", "top_py_lab", "top_pz_lab"],
                library="np",
            )
            delphes_vectors = delphes_tree.arrays(
                [
                    "LargeRJet/LargeRJet.PT",
                    "LargeRJet/LargeRJet.Eta",
                    "LargeRJet/LargeRJet.Phi",
                ],
                library="ak",
            )

            if labels_tree.num_entries != delphes_tree.num_entries:
                raise ValueError(
                    f"Entry mismatch for run {run_index}: labels={labels_tree.num_entries}, "
                    f"delphes={delphes_tree.num_entries}"
                )

            top_eta = _eta_from_momentum(
                label_vectors["top_px_lab"],
                label_vectors["top_py_lab"],
                label_vectors["top_pz_lab"],
            )
            top_phi = np.arctan2(label_vectors["top_py_lab"], label_vectors["top_px_lab"])

            order = ak.argsort(delphes_vectors["LargeRJet/LargeRJet.PT"], axis=1, ascending=False)
            leading_eta = ak.firsts(delphes_vectors["LargeRJet/LargeRJet.Eta"][order])
            leading_phi = ak.firsts(delphes_vectors["LargeRJet/LargeRJet.Phi"][order])

            leading_eta = ak.to_numpy(ak.fill_none(leading_eta, np.nan))
            leading_phi = ak.to_numpy(ak.fill_none(leading_phi, np.nan))

            delta_eta = top_eta - leading_eta
            delta_phi = _delta_phi(top_phi, leading_phi)
            yield np.sqrt(delta_eta**2 + delta_phi**2)


def summarize_alignment(
    labels_dir: Path,
    delphes_dir: Path,
    jet_radius: float | None = None,
) -> AlignmentSummary:
    if jet_radius is None:
        jet_radius = _infer_large_r_radius(delphes_dir)

    labels_dir = labels_dir.resolve()
    delphes_dir = delphes_dir.resolve()

    labels_files = _indexed_root_files(labels_dir, prefix="labels_")
    delphes_files = _indexed_root_files(delphes_dir, prefix="dataset_")

    delta_r = np.concatenate(list(_iter_alignment_delta_r(labels_files, delphes_files)))
    has_large_r_jet = np.isfinite(delta_r)
    matched = has_large_r_jet & (delta_r <= jet_radius)

    total_events = int(delta_r.size)
    events_with_large_r_jet = int(np.count_nonzero(has_large_r_jet))
    events_without_large_r_jet = total_events - events_with_large_r_jet
    matched_events = int(np.count_nonzero(matched))
    unmatched_events_with_large_r_jet = events_with_large_r_jet - matched_events

    if total_events == 0:
        raise ValueError(f"No events found across {labels_dir} and {delphes_dir}")
    if events_with_large_r_jet == 0:
        raise ValueError("No events contain a Large-R jet; cannot evaluate alignment")

    inclusive_match_fraction = matched_events / total_events
    conditional_match_fraction = matched_events / events_with_large_r_jet
    conditional_miss_fraction = unmatched_events_with_large_r_jet / events_with_large_r_jet
    finite_delta_r = delta_r[has_large_r_jet]

    return AlignmentSummary(
        labels_dir=labels_dir,
        delphes_dir=delphes_dir,
        jet_radius=float(jet_radius),
        file_count=len(labels_files),
        total_events=total_events,
        events_with_large_r_jet=events_with_large_r_jet,
        events_without_large_r_jet=events_without_large_r_jet,
        matched_events=matched_events,
        unmatched_events_with_large_r_jet=unmatched_events_with_large_r_jet,
        inclusive_match_fraction=inclusive_match_fraction,
        conditional_match_fraction=conditional_match_fraction,
        conditional_miss_fraction=conditional_miss_fraction,
        max_delta_r_with_large_r_jet=float(np.max(finite_delta_r)),
        p95_delta_r_with_large_r_jet=float(np.quantile(finite_delta_r, 0.95)),
    )


def assert_alignment_tolerances(
    summary: AlignmentSummary,
    *,
    min_match_fraction: float,
    max_no_jet_fraction: float | None,
) -> None:
    if not 0.0 <= min_match_fraction <= 1.0:
        raise ValueError(f"min_match_fraction must be between 0 and 1, got {min_match_fraction}")
    if max_no_jet_fraction is not None and not 0.0 <= max_no_jet_fraction <= 1.0:
        raise ValueError(f"max_no_jet_fraction must be between 0 and 1, got {max_no_jet_fraction}")

    if summary.conditional_match_fraction < min_match_fraction:
        raise AssertionError(
            "Hadronic-top / leading-Large-R-jet alignment fell below tolerance: "
            f"match fraction among events with a Large-R jet is "
            f"{summary.conditional_match_fraction:.4f}, below the required {min_match_fraction:.4f}. "
            f"Matched {summary.matched_events} / {summary.events_with_large_r_jet} events "
            f"using ΔR <= {summary.jet_radius:.4f}."
        )

    if max_no_jet_fraction is not None:
        no_jet_fraction = summary.events_without_large_r_jet / summary.total_events
        if no_jet_fraction > max_no_jet_fraction:
            raise AssertionError(
                "Too many events lack a Large-R jet: "
                f"no-jet fraction is {no_jet_fraction:.4f}, above the allowed {max_no_jet_fraction:.4f}."
            )


def _format_summary(summary: AlignmentSummary) -> str:
    no_jet_fraction = summary.events_without_large_r_jet / summary.total_events
    return "\n".join(
        [
            f"labels dir: {summary.labels_dir}",
            f"delphes dir: {summary.delphes_dir}",
            f"files checked: {summary.file_count}",
            f"Large-R radius used for match test: {summary.jet_radius:.4f}",
            f"total events: {summary.total_events}",
            f"events with >=1 Large-R jet: {summary.events_with_large_r_jet}",
            f"events without a Large-R jet: {summary.events_without_large_r_jet} ({no_jet_fraction:.4%})",
            f"matched events with ΔR(top, leading Large-R jet) <= R: {summary.matched_events}",
            (
                "match fraction among events with a Large-R jet: "
                f"{summary.conditional_match_fraction:.4%}"
            ),
            f"inclusive match fraction over all events: {summary.inclusive_match_fraction:.4%}",
            f"95th percentile ΔR when a Large-R jet exists: {summary.p95_delta_r_with_large_r_jet:.4f}",
            f"max ΔR when a Large-R jet exists: {summary.max_delta_r_with_large_r_jet:.4f}",
        ]
    )


def _cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate that the hadronic-top truth direction falls inside the leading Large-R jet "
            "for most events that contain a Large-R jet."
        )
    )
    parser.add_argument("--labels-dir", required=True, type=Path)
    parser.add_argument("--delphes-dir", required=True, type=Path)
    parser.add_argument(
        "--jet-radius",
        type=float,
        default=None,
        help=(
            "Large-R jet radius to use for the ΔR cone match. "
            "If omitted, the script tries to read LargeRJetRadius from the rendered Delphes card."
        ),
    )
    parser.add_argument(
        "--min-match-fraction",
        type=float,
        default=0.90,
        help=(
            "Minimum allowed fraction of events with a Large-R jet whose hadronic-top truth direction "
            "falls within ΔR <= R of the leading Large-R jet."
        ),
    )
    parser.add_argument(
        "--max-no-jet-fraction",
        type=float,
        default=None,
        help="Optional upper bound on the fraction of events that have no Large-R jet.",
    )
    return parser


def main() -> int:
    args = _cli_parser().parse_args()
    summary = summarize_alignment(
        labels_dir=args.labels_dir,
        delphes_dir=args.delphes_dir,
        jet_radius=args.jet_radius,
    )
    print(_format_summary(summary))
    assert_alignment_tolerances(
        summary,
        min_match_fraction=args.min_match_fraction,
        max_no_jet_fraction=args.max_no_jet_fraction,
    )
    print("Alignment test passed.")
    return 0


def _env_path(name: str) -> Path | None:
    value = os.getenv(name)
    return Path(value) if value else None


def _debug_config_output_dirs() -> tuple[Path, Path]:
    repo_root = Path(__file__).resolve().parents[1]
    config_path = Path(os.getenv("TP_CONFIG", repo_root / "debug.config")).resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Test config not found: {config_path}")

    # Match run_job.sh by sourcing the selected shell config in an isolated
    # subprocess, so computed tags (for example PY_tag from R/minJetpT) stay
    # synchronized with debug.config. The config is project-owned and is
    # already sourced by the production runner.
    result = subprocess.run(
        [
            "bash",
            "-c",
            'set -e; source "$1"; printf "%s\\n" "$MG_tag" "$DELPHES_tag"',
            "alignment-test-config",
            str(config_path),
        ],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    values = result.stdout.splitlines()
    if len(values) != 2 or not all(values):
        raise ValueError(
            f"Expected MG_tag and DELPHES_tag in {config_path}; got output {result.stdout!r}"
        )

    mg_tag, delphes_tag = values
    return (
        repo_root / "madgraph" / f"pp_tt_semi_full_{mg_tag}",
        repo_root / "delphes" / f"WS_{delphes_tag}",
    )


def test_hadronic_top_is_mostly_captured_by_leading_large_r_jet() -> None:
    default_labels_dir, default_delphes_dir = _debug_config_output_dirs()
    labels_dir = _env_path("TP_LABELS_DIR") or default_labels_dir
    delphes_dir = _env_path("TP_DELPHES_DIR") or default_delphes_dir
    jet_radius = os.getenv("TP_LARGE_R_RADIUS")
    min_match_fraction = float(os.getenv("TP_MIN_MATCH_FRACTION", "0.90"))
    max_no_jet_fraction_env = os.getenv("TP_MAX_NO_JET_FRACTION")
    max_no_jet_fraction = (
        float(max_no_jet_fraction_env) if max_no_jet_fraction_env is not None else None
    )

    summary = summarize_alignment(
        labels_dir=labels_dir,
        delphes_dir=delphes_dir,
        jet_radius=float(jet_radius) if jet_radius is not None else None,
    )
    print(_format_summary(summary))
    assert_alignment_tolerances(
        summary,
        min_match_fraction=min_match_fraction,
        max_no_jet_fraction=max_no_jet_fraction,
    )


if __name__ == "__main__":
    raise SystemExit(main())
