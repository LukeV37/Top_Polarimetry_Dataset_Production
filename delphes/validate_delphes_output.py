#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import uproot


def count_hepmc_events(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if len(line) >= 2 and line[0] == "E" and line[1].isspace():
                count += 1
    return count


def flatten_single_value_objects(values: np.ndarray) -> np.ndarray:
    flattened = []
    for value in values:
        if isinstance(value, np.ndarray):
            if value.size != 1:
                raise SystemExit(f"Expected one Event.Number value per entry, got {value!r}")
            flattened.append(int(value.reshape(-1)[0]))
        else:
            flattened.append(int(value))
    return np.asarray(flattened, dtype=np.int64)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Delphes output against its HepMC input.")
    parser.add_argument("--hepmc", required=True, help="Input HepMC file used for Delphes.")
    parser.add_argument("--root", required=True, help="Delphes ROOT output file to validate.")
    args = parser.parse_args()

    hepmc_path = Path(args.hepmc)
    root_path = Path(args.root)

    if not hepmc_path.is_file() or hepmc_path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty HepMC input: {hepmc_path}")
    if not root_path.is_file() or root_path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty Delphes ROOT output: {root_path}")

    hepmc_events = count_hepmc_events(hepmc_path)

    with uproot.open(root_path) as root_file:
        if "Delphes" not in root_file:
            raise SystemExit(f"Delphes tree not found in {root_path}")

        tree = root_file["Delphes"]
        event_numbers = flatten_single_value_objects(tree["Event/Event.Number"].array(library="np"))

    root_events = len(event_numbers)

    if root_events != hepmc_events:
        raise SystemExit(
            f"Delphes/HepMC event-count mismatch for {root_path}: "
            f"Delphes has {root_events}, HepMC has {hepmc_events}"
        )

    expected_sequence = np.arange(root_events, dtype=event_numbers.dtype)
    if not np.array_equal(event_numbers, expected_sequence):
        preview = ", ".join(str(int(x)) for x in event_numbers[:10])
        raise SystemExit(
            f"Non-sequential Event.Number values in {root_path}. "
            f"Expected 0..{root_events - 1}; first values: [{preview}]"
        )

    print(f"Validated {root_events} Delphes events against {hepmc_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
