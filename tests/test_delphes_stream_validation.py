"""Exercise stream validation without running Pythia or Delphes."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import uproot


VALIDATOR = Path(__file__).resolve().parents[1] / "delphes/validate_delphes_output.py"


@pytest.mark.parametrize(
    ("completion", "numbers", "error"),
    [
        ("HEPMC_EVENTS_WRITTEN=2", [0, 1], None),
        ("Successfully wrote 2 events to /proc/self/fd/3", [0, 1], None),
        ("", [0, 1], "Expected one successful, nonempty HepMC stream"),
        ("HEPMC_EVENTS_WRITTEN=0", [], "Expected one successful, nonempty HepMC stream"),
        ("HEPMC_EVENTS_WRITTEN=1", [0], "MadGraph/HepMC event-count mismatch"),
        ("HEPMC_EVENTS_WRITTEN=2", [0], "Delphes/HepMC event-count mismatch"),
        ("HEPMC_EVENTS_WRITTEN=2", [1, 0], "Non-sequential Event.Number"),
    ],
)
def test_stream_validation(tmp_path, completion, numbers, error):
    log = tmp_path / "pythia.log"
    log.write_text(f"Pythia diagnostics\n{completion}\n")
    labels = tmp_path / "labels.root"
    with uproot.recreate(labels) as root:
        root.mktree("labels", {"label": np.arange(2)})
    reco = tmp_path / "delphes.root"
    with uproot.recreate(reco) as root:
        root.mktree("Delphes", {"Event/Event.Number": np.asarray(numbers, dtype=np.int64)})

    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "--pythia-log", str(log), "--labels", str(labels), "--root", str(reco)],
        capture_output=True, text=True,
    )
    if error is None:
        assert result.returncode == 0, result.stderr
        assert "Validated 2 Delphes events" in result.stdout
    else:
        assert result.returncode != 0
        assert error in result.stderr
