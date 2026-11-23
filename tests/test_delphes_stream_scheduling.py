"""Test streamed-run scheduling with small stand-ins, without generating events."""

import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest


REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def fake_pipeline(tmp_path):
    delphes = tmp_path / "delphes"
    pythia = tmp_path / "pythia/src"
    binary_dir = tmp_path / "submodules/delphes-v3.5.0"
    labels = tmp_path / "madgraph/pp_tt_semi_full_test"
    for path in (delphes, pythia, binary_dir, labels):
        path.mkdir(parents=True)

    (delphes / "run.sh").symlink_to(REPO / "delphes/run.sh")
    (delphes / "render_card.py").write_text(
        "import pathlib, sys\n"
        "pathlib.Path(sys.argv[sys.argv.index('--output') + 1]).write_text('card')\n"
    )
    (delphes / "validate_delphes_output.py").write_text(
        "import os, sys\n"
        "sys.exit(1 if os.environ.get('FAIL_VALIDATION') else 0)\n"
    )
    (pythia / "makefile").write_text("generate_dataset:\n\t@true\n")
    producer = pythia / "run_dataset"
    producer.write_text(
        "#!/usr/bin/env python3\n"
        "import os, pathlib, sys, time\n"
        "run = sys.argv[3]\n"
        "gate = pathlib.Path(os.environ['GATE'])\n"
        "(gate / ('started_' + run)).touch()\n"
        "while (gate / ('block_' + run)).exists(): time.sleep(0.01)\n"
        "if os.environ.get('FAIL_PRODUCER') == run: sys.exit(1)\n"
        "with open(sys.argv[4], 'w') as stream: stream.write('event\\n')\n"
        "print('HEPMC_EVENTS_WRITTEN=1')\n"
    )
    producer.chmod(0o755)
    consumer = binary_dir / "DelphesHepMC2"
    consumer.write_text(
        "#!/usr/bin/env python3\n"
        "import os, pathlib, sys\n"
        "sys.stdin.read() if sys.argv[3] == '-' else pathlib.Path(sys.argv[3]).read_text()\n"
        "if os.environ.get('FAIL_CONSUMER'): sys.exit(1)\n"
        "pathlib.Path(sys.argv[2]).write_text('root')\n"
    )
    consumer.chmod(0o755)
    for run in range(10):
        (labels / f"labels_test_{run}.root").write_text("label")

    gate = tmp_path / "gate"
    gate.mkdir()
    env = dict(os.environ, GATE=str(gate), PATH=f"{Path(sys.executable).parent}:{os.environ['PATH']}")
    return delphes, gate, env


def wait_for_files(paths, timeout=15):
    deadline = time.monotonic() + timeout
    while not all(path.exists() for path in paths):
        if time.monotonic() >= deadline:
            raise AssertionError(f"Timed out waiting for: {[str(p) for p in paths if not p.exists()]}")
        time.sleep(0.01)


def run_pipeline(delphes, env, runs, workers, stream=True):
    return subprocess.Popen(
        ["bash", "run.sh", "test", "test", str(runs), str(workers)] + (["test"] if stream else []),
        cwd=delphes, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, start_new_session=True,
    )


def stop_pipeline(process):
    if process.poll() is None:
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()


def test_ten_streamed_runs_start_together(fake_pipeline):
    delphes, gate, env = fake_pipeline
    for run in range(10):
        (gate / f"block_{run}").touch()
    process = run_pipeline(delphes, env, 10, 10)
    try:
        wait_for_files([gate / f"started_{run}" for run in range(10)])
        for run in range(10):
            (gate / f"block_{run}").unlink()
        stdout, stderr = process.communicate(timeout=15)
        assert process.returncode == 0, (stdout, stderr)
        assert all((delphes / f"WS_test/dataset_test_{run}.root").is_file() for run in range(10))
    finally:
        stop_pipeline(process)


def test_rejects_more_runs_than_workers(fake_pipeline):
    delphes, gate, env = fake_pipeline
    process = run_pipeline(delphes, env, 3, 2)
    try:
        stdout, stderr = process.communicate(timeout=15)
        assert process.returncode != 0
        assert "num_runs <= max_cpu_cores" in stderr
        assert not list(gate.glob("started_*"))
        assert not (delphes / "WS_test").exists()
    finally:
        stop_pipeline(process)


def test_top_level_rejects_too_many_runs_before_generation(tmp_path):
    config = tmp_path / "job.config"
    config.write_text(
        "stream_hepmc=true\nbypass_madgraph=false\nbypass_pythia=false\n"
        "bypass_delphes=false\nnum_runs=11\nmax_cpu_cores=10\n"
    )
    result = subprocess.run(
        ["bash", str(REPO / "run_job.sh"), str(config)], cwd=tmp_path,
        capture_output=True, text=True, timeout=15,
    )
    assert result.returncode != 0
    assert "num_runs <= max_cpu_cores" in result.stderr


def test_file_based_runs_still_use_batches(fake_pipeline):
    delphes, _, env = fake_pipeline
    inputs = delphes.parent / "pythia/WS_test"
    inputs.mkdir()
    for run in range(3):
        (inputs / f"dataset_test_{run}.hepmc").write_text("event")
    process = run_pipeline(delphes, env, 3, 2, stream=False)
    try:
        stdout, stderr = process.communicate(timeout=15)
        assert process.returncode == 0, (stdout, stderr)
        assert "batch 1" in stdout
        assert all((delphes / f"WS_test/dataset_test_{run}.root").is_file() for run in range(3))
    finally:
        stop_pipeline(process)


def test_validation_failure_does_not_publish_output(fake_pipeline):
    delphes, _, env = fake_pipeline
    env["FAIL_VALIDATION"] = "1"
    process = run_pipeline(delphes, env, 1, 1)
    try:
        stdout, stderr = process.communicate(timeout=15)
        assert process.returncode != 0, (stdout, stderr)
        assert not (delphes / "WS_test/dataset_test_0.root").exists()
        assert not list((delphes / "WS_test").glob(".stream_*"))
    finally:
        stop_pipeline(process)


@pytest.mark.parametrize("failure", ["FAIL_PRODUCER", "FAIL_CONSUMER"])
def test_pipeline_failure_does_not_publish_output(fake_pipeline, failure):
    delphes, gate, env = fake_pipeline
    env[failure] = "0" if failure == "FAIL_PRODUCER" else "1"
    process = run_pipeline(delphes, env, 1, 1)
    try:
        stdout, stderr = process.communicate(timeout=15)
        assert process.returncode != 0, (stdout, stderr)
        assert not (delphes / "WS_test/dataset_test_0.root").exists()
        assert not list((delphes / "WS_test").glob(".stream_*"))
        assert not list((delphes / "WS_test").glob(".queue.*"))
    finally:
        stop_pipeline(process)
