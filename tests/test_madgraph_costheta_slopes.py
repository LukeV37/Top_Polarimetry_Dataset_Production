"""Regression checks for cos(theta) slopes in the MadGraph truth labels.

This test generates small MadGraph samples for left-, right-, and unpolarized
tops, runs the repository's existing LHE-to-label pipeline, and checks the
resulting ``costheta_down`` and ``costheta_bottom`` distributions.

The fitted quantity is the slope parameter ``beta`` in the normalized linear
shape

    f(cos(theta)) = 0.5 + 0.5 * beta * cos(theta)

which is obtained by fitting a straight line to the normalized histogram over
``[-1, 1]`` and converting the fitted line slope ``m`` to ``beta = 2m``.

Expected values:

- U: down = 0.0, bottom = 0.0
- L: down = -1.0, bottom = +0.4
- R: down = +1.0, bottom = -0.4

Because this is an integration test that runs MadGraph, it is intentionally
configurable through environment variables for local use:

- ``TOP_POLARIMETRY_COSTHETA_TEST_EVENTS`` (default: ``10000``)
- ``TOP_POLARIMETRY_COSTHETA_TEST_TOP_PT_CUT`` (default: ``1``)
- ``TOP_POLARIMETRY_COSTHETA_TEST_BINS`` (default: ``60``)
- ``TOP_POLARIMETRY_COSTHETA_TEST_SLOPE_TOL`` (default: ``0.10``)
- ``TOP_POLARIMETRY_COSTHETA_TEST_TIMEOUT`` in seconds (default: ``1200``)
- ``TOP_POLARIMETRY_KEEP_TEST_ARTIFACTS=1`` to keep generated outputs
"""

from __future__ import annotations

import os
import shutil
import subprocess
import textwrap
import uuid
from pathlib import Path
from typing import Iterator

import numpy as np
import pytest
import uproot


REPO_ROOT = Path(__file__).resolve().parents[1]
MADGRAPH_DIR = REPO_ROOT / "madgraph"

NUM_EVENTS = int(os.environ.get("TOP_POLARIMETRY_COSTHETA_TEST_EVENTS", "10000"))
TOP_PT_CUT = float(os.environ.get("TOP_POLARIMETRY_COSTHETA_TEST_TOP_PT_CUT", "1"))
FIT_BINS = int(os.environ.get("TOP_POLARIMETRY_COSTHETA_TEST_BINS", "60"))
SLOPE_TOLERANCE = float(
    os.environ.get("TOP_POLARIMETRY_COSTHETA_TEST_SLOPE_TOL", "0.10")
)
INTERCEPT_TOLERANCE = float(
    os.environ.get("TOP_POLARIMETRY_COSTHETA_TEST_INTERCEPT_TOL", "0.08")
)
TIMEOUT_SECONDS = int(
    os.environ.get("TOP_POLARIMETRY_COSTHETA_TEST_TIMEOUT", "1200")
)
KEEP_ARTIFACTS = os.environ.get("TOP_POLARIMETRY_KEEP_TEST_ARTIFACTS") == "1"

EXPECTED_SLOPES = {
    "U": {"costheta_down": 0.0, "costheta_bottom": 0.0},
    "L": {"costheta_down": -1.0, "costheta_bottom": 0.4},
    "R": {"costheta_down": 1.0, "costheta_bottom": -0.4},
}


def _fit_histogram_beta(values: np.ndarray, bins: int = FIT_BINS) -> tuple[float, float]:
    """Fit a line to the normalized histogram and return (beta, intercept)."""
    values = np.asarray(values, dtype=np.float64)
    assert values.ndim == 1
    assert len(values) > 0, "Cannot fit an empty cos(theta) sample"

    counts, edges = np.histogram(values, bins=bins, range=(-1.0, 1.0))
    assert counts.sum() == len(values), "Histogram lost entries outside [-1, 1]"

    bin_width = edges[1] - edges[0]
    density = counts / (counts.sum() * bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Weight by approximate inverse Poisson uncertainty. Empty bins contribute
    # no information to the weighted fit.
    weights = np.sqrt(counts.astype(np.float64))
    nonzero = weights > 0
    assert nonzero.sum() >= 2, "Need at least two populated bins for a line fit"

    slope, intercept = np.polyfit(
        centers[nonzero], density[nonzero], deg=1, w=weights[nonzero]
    )
    beta = 2.0 * slope
    return float(beta), float(intercept)


@pytest.fixture(params=[pytest.param("L"), pytest.param("R"), pytest.param("U")])
def generated_sample(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Iterator[tuple[str, str, Path]]:
    polarization = str(request.param)
    tag = f"pytest_costheta_{polarization.lower()}_{uuid.uuid4().hex[:10]}"

    job_config = tmp_path / f"job_{polarization}.config"
    job_config.write_text(
        textwrap.dedent(
            f"""
            bypass_madgraph=false
            bypass_pythia=true
            bypass_preprocessing=true
            bypass_train=true
            process=\"{polarization}\"
            generation=\"inclusive\"
            MG_tag=\"{tag}\"
            num_runs=1
            num_events_per_run=\"{NUM_EVENTS}\"
            max_cpu_cores=1
            seed=12
            top_pT_cut={TOP_PT_CUT:g}
            """
        ).strip()
        + "\n"
    )

    result = subprocess.run(
        [str(REPO_ROOT / "run_job.sh"), str(job_config)],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=TIMEOUT_SECONDS,
        check=False,
    )
    if result.returncode != 0:
        pytest.fail(
            f"MadGraph job failed for polarization {polarization} "
            f"(exit {result.returncode}). "
            "Generated files were left in place for diagnosis.\n"
            f"--- job output (last 30000 characters) ---\n{result.stdout[-30000:]}"
        )

    process_dir = MADGRAPH_DIR / f"pp_tt_semi_full_{tag}"
    yield polarization, tag, process_dir

    if not KEEP_ARTIFACTS:
        shutil.rmtree(process_dir, ignore_errors=True)
        (MADGRAPH_DIR / f"MadGraph_{tag}.log").unlink(missing_ok=True)


def test_costheta_slopes_match_expected_polarization(
    generated_sample: tuple[str, str, Path],
) -> None:
    polarization, tag, process_dir = generated_sample
    labels_root = process_dir / f"labels_{tag}_0.root"
    assert labels_root.is_file(), f"Missing labels ROOT file: {labels_root}"

    with uproot.open(labels_root) as label_file:
        labels = label_file["labels"]
        arrays = labels.arrays(["costheta_down", "costheta_bottom"], library="np")

    assert len(arrays["costheta_down"]) == len(arrays["costheta_bottom"]) > 0

    for branch, expected_beta in EXPECTED_SLOPES[polarization].items():
        values = arrays[branch]
        assert np.isfinite(values).all(), f"Non-finite values in {branch} for {labels_root}"
        assert ((values >= -1.00001) & (values <= 1.00001)).all(), (
            f"{branch} contains values outside [-1, 1] in {labels_root}"
        )

        fitted_beta, intercept = _fit_histogram_beta(values)
        print(
            f"{polarization} {branch}: fitted beta={fitted_beta:+.4f} "
            f"(expected {expected_beta:+.4f} ± {SLOPE_TOLERANCE:.3f}), "
            f"histogram line slope={fitted_beta / 2.0:+.4f}, "
            f"intercept={intercept:.4f}, entries={len(values)}"
        )

        assert intercept == pytest.approx(0.5, abs=INTERCEPT_TOLERANCE), (
            f"Unexpected normalized-histogram intercept for {polarization} {branch}: "
            f"expected 0.5 ± {INTERCEPT_TOLERANCE}, got {intercept:.4f}"
        )
        assert fitted_beta == pytest.approx(expected_beta, abs=SLOPE_TOLERANCE), (
            f"Unexpected cos(theta) slope for {polarization} {branch}: "
            f"expected {expected_beta:+.3f} ± {SLOPE_TOLERANCE}, "
            f"got {fitted_beta:+.4f}"
        )
