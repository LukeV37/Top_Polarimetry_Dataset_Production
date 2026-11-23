#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


ALGORITHM_CODES = {
    "kt": (4, "kt"),
    "cambridge": (5, "cambridge"),
    "cambridge_aachen": (5, "cambridge"),
    "ca": (5, "cambridge"),
    "antikt": (6, "antikt"),
}


def normalize_algorithm(name: str) -> tuple[int, str]:
    key = name.strip().lower().replace("-", "").replace("/", "_").replace(" ", "_")
    key = key.replace("antikt", "antikt").replace("antik_t", "antikt")
    key = key.replace("cambridgeaachen", "cambridge_aachen")
    if key == "anti_kt":
        key = "antikt"
    if key == "cambridge/aachen":
        key = "cambridge_aachen"

    if key not in ALGORITHM_CODES:
        choices = ", ".join(sorted(ALGORITHM_CODES))
        raise ValueError(f"Unsupported jet algorithm '{name}'. Supported values: {choices}")

    return ALGORITHM_CODES[key]


def positive_float(value: str, *, allow_zero: bool = False) -> float:
    parsed = float(value)
    if allow_zero:
        if parsed < 0.0:
            raise ValueError(f"Expected a non-negative value, got {value}")
    elif parsed <= 0.0:
        raise ValueError(f"Expected a positive value, got {value}")
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description="Render the Top Polarimetry Delphes ATLAS card.")
    parser.add_argument("--template", required=True, help="Path to the Delphes card template.")
    parser.add_argument("--output", required=True, help="Path to write the rendered Delphes card.")
    parser.add_argument("--small-r-algorithm", required=True)
    parser.add_argument("--small-r", required=True)
    parser.add_argument("--small-r-ptmin", required=True)
    parser.add_argument("--large-r-algorithm", required=True)
    parser.add_argument("--large-r", required=True)
    parser.add_argument("--large-r-ptmin", required=True)
    parser.add_argument("--small-r-jet-flavor-dr", required=True)
    parser.add_argument("--small-r-parton-ptmin", required=True)
    parser.add_argument("--small-r-parton-eta-max", required=True)
    args = parser.parse_args()

    small_r_algorithm_code, small_r_algorithm_name = normalize_algorithm(args.small_r_algorithm)
    large_r_algorithm_code, large_r_algorithm_name = normalize_algorithm(args.large_r_algorithm)

    small_r = positive_float(args.small_r)
    small_r_ptmin = positive_float(args.small_r_ptmin, allow_zero=True)
    large_r = positive_float(args.large_r)
    large_r_ptmin = positive_float(args.large_r_ptmin, allow_zero=True)
    small_r_jet_flavor_dr = positive_float(args.small_r_jet_flavor_dr)
    small_r_parton_ptmin = positive_float(args.small_r_parton_ptmin, allow_zero=True)
    small_r_parton_eta_max = positive_float(args.small_r_parton_eta_max)

    template_path = Path(args.template)
    output_path = Path(args.output)

    rendered = template_path.read_text()
    replacements = {
        "__TP_SMALL_R_JET_ALGO_NAME__": small_r_algorithm_name,
        "__TP_SMALL_R_JET_ALGO_CODE__": str(small_r_algorithm_code),
        "__TP_SMALL_R__": f"{small_r}",
        "__TP_SMALL_R_JET_PTMIN__": f"{small_r_ptmin}",
        "__TP_LARGE_R_JET_ALGO_NAME__": large_r_algorithm_name,
        "__TP_LARGE_R_JET_ALGO_CODE__": str(large_r_algorithm_code),
        "__TP_LARGE_R__": f"{large_r}",
        "__TP_LARGE_R_JET_PTMIN__": f"{large_r_ptmin}",
        "__TP_SMALL_R_JET_FLAVOR_DR__": f"{small_r_jet_flavor_dr}",
        "__TP_SMALL_R_PARTON_PTMIN__": f"{small_r_parton_ptmin}",
        "__TP_SMALL_R_PARTON_ETA_MAX__": f"{small_r_parton_eta_max}",
    }

    for placeholder, value in replacements.items():
        rendered = rendered.replace(placeholder, value)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
