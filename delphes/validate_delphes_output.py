#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import re

import numpy as np
import uproot


def count_hepmc_events(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if len(line) >= 2 and line[0] == "E" and line[1].isspace():
                count += 1
    return count


def count_streamed_events(log: Path) -> int:
    # Pythia emits this completion line only after closing the HepMC stream.
    text = log.read_text(encoding="utf-8", errors="replace")
    counts = re.findall(r"^HEPMC_EVENTS_WRITTEN=(\d+)$", text, re.MULTILINE)
    if not counts:  # A running job may still use the previous Pythia executable.
        counts = re.findall(r"^Successfully wrote (\d+) events to .+$", text, re.MULTILINE)
    if len(counts) != 1 or int(counts[0]) == 0:
        raise SystemExit(f"Expected one successful, nonempty HepMC stream in {log}")
    return int(counts[0])


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
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--hepmc", type=Path, help="Input HepMC file used for Delphes.")
    source.add_argument("--pythia-log", type=Path, help="Pythia log for a streamed HepMC input.")
    parser.add_argument("--labels", type=Path, help="MadGraph label ROOT file (required for streaming).")
    parser.add_argument("--root", required=True, help="Delphes ROOT output file to validate.")
    args = parser.parse_args()
    root_path = Path(args.root)

    if args.pythia_log and not args.labels:
        parser.error("--pythia-log requires --labels")
    if not root_path.is_file() or root_path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty Delphes ROOT output: {root_path}")

    if args.hepmc:
        if not args.hepmc.is_file() or args.hepmc.stat().st_size == 0:
            raise SystemExit(f"Missing or empty HepMC input: {args.hepmc}")
        hepmc_events = count_hepmc_events(args.hepmc)
    else:
        hepmc_events = count_streamed_events(args.pythia_log)

    if args.labels:
        with uproot.open(args.labels) as label_file:
            label_events = label_file["labels"].num_entries
        if label_events != hepmc_events:
            raise SystemExit(f"MadGraph/HepMC event-count mismatch: labels={label_events}, HepMC={hepmc_events}")

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

    print(f"Validated {root_events} Delphes events against {(args.hepmc or args.pythia_log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
