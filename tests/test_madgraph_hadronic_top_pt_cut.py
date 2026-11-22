"""Regression checks for the MadGraph hadronic-top pT cut.

The integration checks create temporary job configurations and run a small
MadGraph job at each threshold. They require the repository's configured
MadGraph installation and Python/ROOT dependencies.
"""

from __future__ import annotations

import gzip
import math
import re
import shutil
import subprocess
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterator

import numpy as np
import pytest
import uproot


REPO_ROOT = Path(__file__).resolve().parents[1]
MADGRAPH_DIR = REPO_ROOT / "madgraph"

SAMPLES = [
    pytest.param(200.0, id="200-GeV"),
    pytest.param(300.0, id="300-GeV"),
]


def _normalise_fortran(source: str) -> str:
    """Normalise case and whitespace for checks against MadGraph output."""
    return re.sub(r"\s+", " ", source.upper()).strip()


def _hadronic_top_pts(lhe_path: Path) -> Iterator[float]:
    """Yield lab-frame pT from the final b and two light-quark daughters."""
    with gzip.open(lhe_path, "rb") as source:
        for _, element in ET.iterparse(source, events=("end",)):
            if element.tag.rsplit("}", 1)[-1] != "event":
                continue

            fields = (element.text or "").split()
            assert fields, f"Empty LHE event in {lhe_path}"
            n_particles = int(fields[0])
            expected_fields = 6 + 13 * n_particles
            assert len(fields) >= expected_fields, (
                f"Truncated LHE event in {lhe_path}: expected at least "
                f"{expected_fields} fields, found {len(fields)}"
            )
            particles = [
                fields[6 + 13 * index : 6 + 13 * (index + 1)]
                for index in range(n_particles)
            ]
            final_state = [particle for particle in particles if int(particle[1]) == 1]
            hadronic_bottom = [
                particle for particle in final_state if int(particle[0]) == 5
            ]
            light_quarks = [
                particle
                for particle in final_state
                if 1 <= abs(int(particle[0])) <= 4
            ]
            assert len(hadronic_bottom) == 1 and len(light_quarks) == 2, (
                f"Expected one hadronic b and two light-quark daughters in "
                f"{lhe_path}; found b={len(hadronic_bottom)}, "
                f"light quarks={len(light_quarks)}"
            )

            daughters = hadronic_bottom + light_quarks
            px = sum(float(particle[6]) for particle in daughters)
            py = sum(float(particle[7]) for particle in daughters)
            yield math.hypot(px, py)
            element.clear()


@pytest.fixture(params=SAMPLES)
def generated_sample(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Iterator[tuple[float, str, Path]]:
    threshold = float(request.param)
    tag = f"pytest_hadtoppt_{int(threshold)}_{uuid.uuid4().hex[:10]}"
    job_config = tmp_path / f"job_{int(threshold)}.config"
    job_config.write_text(
        "\n".join(
            (
                "bypass_madgraph=false",
                "bypass_pythia=true",
                "bypass_preprocessing=true",
                "bypass_train=true",
                'process="U"',
                'generation="inclusive"',
                f'MG_tag="{tag}"',
                "num_runs=1",
                'num_events_per_run="100"',
                "max_cpu_cores=1",
                "seed=12",
                f"top_pT_cut={threshold:g}",
                "",
            )
        )
    )

    result = subprocess.run(
        [str(REPO_ROOT / "run_job.sh"), str(job_config)],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=600,
        check=False,
    )
    if result.returncode != 0:
        pytest.fail(
            f"MadGraph job failed for {threshold:g} GeV (exit {result.returncode}). "
            "The generated files and log were left in place for diagnosis.\n"
            f"--- job output (last 30000 characters) ---\n{result.stdout[-30000:]}"
        )

    process_dir = MADGRAPH_DIR / f"pp_tt_semi_full_{tag}"
    yield threshold, tag, process_dir

    # Keep artifacts when assertions fail; remove only this test's uniquely
    # tagged outputs after a successful validation.
    shutil.rmtree(process_dir)
    log_file = MADGRAPH_DIR / f"MadGraph_{tag}.log"
    log_file.unlink(missing_ok=True)


def test_custom_cut_source_uses_configured_hadronic_top_pt() -> None:
    source_path = MADGRAPH_DIR / "configs" / "hadronic_top_pt_cut.f"
    source = _normalise_fortran(source_path.read_text())

    assert "LOGICAL FUNCTION DUMMY_CUTS(P)" in source
    assert "PX_HAD_TOP = P(1,3) + P(1,4) + P(1,5)" in source
    assert "PY_HAD_TOP = P(2,3) + P(2,4) + P(2,5)" in source
    assert "PT_HAD_TOP = DSQRT(PX_HAD_TOP**2 + PY_HAD_TOP**2)" in source
    assert "IF (PT_HAD_TOP .LT. MIN_HADRONIC_TOP_PT) THEN" in source
    assert "DUMMY_CUTS = .FALSE." in source


def test_generated_sample_applies_cut_and_matches_truth_labels(
    generated_sample: tuple[float, str, Path],
) -> None:
    threshold, tag, process_dir = generated_sample
    events_dir = process_dir / "Events"

    run_banner = events_dir / "run_01_banner.txt"
    assert run_banner.is_file(), f"Missing MadGraph banner: {run_banner}"
    banner = run_banner.read_text(errors="replace")

    threshold_lines = re.findall(
        r"^\s*([0-9]+(?:\.[0-9]+)?)\s*=\s*min_hadronic_top_pt\b",
        banner,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    assert threshold_lines, f"No min_hadronic_top_pt entry in {run_banner}"
    assert float(threshold_lines[-1]) == pytest.approx(threshold)

    hook_lines = re.findall(
        r"^\s*(.*?)\s*=\s*custom_fcts\b",
        banner,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    assert hook_lines, f"No custom_fcts entry in {run_banner}"
    assert Path(hook_lines[-1].strip()).name == "hadronic_top_pt_cut.f"

    dummy_fct = process_dir / "SubProcesses" / "dummy_fct.f"
    assert dummy_fct.is_file(), f"Missing generated custom function: {dummy_fct}"
    generated_function = _normalise_fortran(dummy_fct.read_text())
    assert "LOGICAL FUNCTION DUMMY_CUTS(P)" in generated_function
    assert "PT_HAD_TOP = DSQRT(PX_HAD_TOP**2 + PY_HAD_TOP**2)" in generated_function
    assert "IF (PT_HAD_TOP .LT. MIN_HADRONIC_TOP_PT) THEN" in generated_function

    run_include = process_dir / "Source" / "run.inc"
    assert run_include.is_file(), f"Missing generated run.inc: {run_include}"
    assert "MIN_HADRONIC_TOP_PT" in run_include.read_text().upper(), (
        "The configured threshold was not declared for the custom Fortran hook"
    )

    generated_cuts = process_dir / "SubProcesses" / "cuts.f"
    assert generated_cuts.is_file(), f"Missing generated cuts.f: {generated_cuts}"
    assert "PT_HAD_TOP" not in generated_cuts.read_text().upper(), (
        "The custom pT cut should be installed through dummy_fct.f, not cuts.f"
    )

    run_indices = sorted(
        int(match.group(1))
        for path in events_dir.iterdir()
        if (match := re.fullmatch(r"run_01_(\d+)", path.name)) and path.is_dir()
    )
    assert run_indices, f"No per-run event directories found in {events_dir}"
    assert run_indices == list(range(len(run_indices))), (
        f"Expected contiguous run indices in {events_dir}; got {run_indices}"
    )

    for run_index in run_indices:
        lhe_file = events_dir / f"run_01_{run_index}" / "unweighted_events.lhe.gz"
        hard_root = process_dir / f"hard_process_{tag}_{run_index}.root"
        labels_root = process_dir / f"labels_{tag}_{run_index}.root"
        assert lhe_file.is_file(), f"Missing LHE file: {lhe_file}"
        assert hard_root.is_file(), f"Missing hard-process ROOT file: {hard_root}"
        assert labels_root.is_file(), f"Missing labels ROOT file: {labels_root}"

        lab_top_pts = np.asarray(list(_hadronic_top_pts(lhe_file)), dtype=np.float64)
        assert len(lab_top_pts) > 0, f"No events found in {lhe_file}"

        with uproot.open(hard_root) as hard_file:
            hard_entries = hard_file["events"].num_entries
        with uproot.open(labels_root) as label_file:
            labels = label_file["labels"]
            label_entries = labels.num_entries
            label_top_pts = labels["top_pT"].array(library="np")
            costheta_down = labels["costheta_down"].array(library="np")
            costheta_bottom = labels["costheta_bottom"].array(library="np")

        assert len(lab_top_pts) == hard_entries == label_entries, (
            f"Entry mismatch for {tag} run {run_index}: "
            f"LHE={len(lab_top_pts)}, hard-process={hard_entries}, "
            f"labels={label_entries}"
        )
        assert np.isfinite(label_top_pts).all(), f"Non-finite top_pT in {labels_root}"
        assert np.isfinite(costheta_down).all() and np.isfinite(costheta_bottom).all(), (
            f"Non-finite angular labels in {labels_root}"
        )
        assert ((costheta_down >= -1.00001) & (costheta_down <= 1.00001)).all()
        assert ((costheta_bottom >= -1.00001) & (costheta_bottom <= 1.00001)).all()

        np.testing.assert_allclose(
            label_top_pts,
            lab_top_pts,
            rtol=1e-6,
            atol=1e-3,
            err_msg=f"Label order or hadronic-top pT mismatch in {tag} run {run_index}",
        )
        assert (lab_top_pts >= threshold - 1e-3).all(), (
            f"LHE has a hadronic-top pT below {threshold:g} GeV in "
            f"{tag} run {run_index}; minimum={lab_top_pts.min():.6f} GeV"
        )
        assert (label_top_pts >= threshold - 1e-3).all(), (
            f"Truth labels have a top_pT below {threshold:g} GeV in "
            f"{tag} run {run_index}; minimum={label_top_pts.min():.6f} GeV"
        )
