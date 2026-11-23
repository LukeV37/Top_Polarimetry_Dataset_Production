# ROOT analysis selection

The analysis reads matching Delphes and MadGraph label files, selects events,
and writes one ROOT file per job. Each output contains one TTree named
`selected`, with reconstructed quantities and the corresponding MadGraph
truth labels. Run the configured pipeline from the repository root:

```bash
./run_job.sh debug.config
```

Set the bypass options in the config to reuse existing inputs. The default
output is
`analysis/WS_<DELPHES_TAG>/dataset_selected_<DELPHES_TAG>_<i>.root`, with a
per-job log under `logs/`. `ANALYSIS_OUTPUT_DIR` can override the output
directory. The runner does not merge the files.

Input files are paired by their `_i` suffix and joined by entry position. The
selection checks equal entry counts and sequential Delphes `Event.Number`
values; it reads both inputs even for events rejected by the physics cuts.
The label tree has no independent event identifier, so this check relies on
upstream stages preserving the same event order. Errors stop the job, and a
failed job does not publish a partial output file.

## Selection overview

Thresholds and capacities are configurable in `job.config` and `debug.config`.
The current defaults are:

1. Select the highest-pT qualifying electron or muon from Delphes' already
   isolated collections: `pT >= 30 GeV` and `|eta| <= 3`. Exact pT ties favor
   the first electron, then the original collection order. No additional
   isolation cut is applied here.
2. Require `MET >= 30 GeV`.
3. Require exactly one large-R jet with `pT >= 250 GeV`, counting jets before
   applying the selected jet's `|eta| <= 3` requirement.
4. Read its charged and neutral energy-flow constituent references without
   additional constituent pT or eta cuts.
5. Discard small-R jets with `DeltaR < 1.5` from the large-R jet axis. This
   threshold is configurable independently of the jet radii; jets are not
   reclustered.
6. Require at least one remaining small-R jet with Delphes b-tag bit 0 set,
   `pT > 30 GeV`, `|eta| < 3`, and `0.3 < DeltaR(lepton, jet) < 3`.

`DeltaR` uses pseudorapidity and the wrapped difference in azimuth. The
Delphes card determines object reconstruction, isolation and tagging; the
above analysis thresholds are not a complete description of those definitions.
Per-job logs record the effective analysis settings, jet radii, cutflow and
aggregate overflow counts. Selection choices remain subject to physics review.

## Output branches

Each tree entry is one accepted event. The tables specify every branch produced
by the current implementation. `Float_t` is a 32-bit float and `Int_t` is a
32-bit signed integer. `Nc` and `Nj` denote the configured fixed array
capacities, 100 and 10 by default. In the ROOT file the array lengths are
concrete numbers, not variable-length arrays.

Momenta, energies and masses are in GeV, azimuthal angles are in radians,
and impact parameters and their errors are in mm. Pseudorapidity, ratios,
counts, IDs and angular cosines are dimensionless.

### Reconstructed lepton, missing energy and large-R jet

| Branch | ROOT type | Units | Description |
|---|---|---|---|
| `lepton_pt` | `Float_t` | GeV | Selected lepton transverse momentum. |
| `lepton_eta` | `Float_t` | — | Selected lepton pseudorapidity. |
| `lepton_phi` | `Float_t` | rad | Selected lepton azimuth. |
| `lepton_charge` | `Int_t` | — | Selected lepton charge in units of the elementary charge. |
| `lepton_flavor` | `Int_t` | — | Selected lepton flavor: `11` for electron or `13` for muon; the sign is not encoded here. |
| `lepton_isolation` | `Float_t` | — | Delphes `IsolationVar`: delta-beta-corrected isolation sum divided by the lepton pT. |
| `lepton_isolation_rho_corr` | `Float_t` | — | Delphes `IsolationVarRhoCorr`: rho-corrected isolation sum divided by the lepton pT. |
| `met` | `Float_t` | GeV | Delphes missing transverse energy magnitude. |
| `met_phi` | `Float_t` | rad | Delphes missing transverse energy azimuth. |
| `large_r_jet_pt` | `Float_t` | GeV | Selected large-R jet transverse momentum. |
| `large_r_jet_eta` | `Float_t` | — | Selected large-R jet pseudorapidity. |
| `large_r_jet_phi` | `Float_t` | rad | Selected large-R jet azimuth. |
| `large_r_jet_mass` | `Float_t` | GeV | Selected large-R jet mass. |

### Large-R jet constituents

These arrays describe the selected jet's Delphes `EFlowTrack`, `EFlowPhoton`
and `EFlowNeutralHadron` references. Entries are sorted by descending stored
`constituent_pt`, preserving the original constituent-reference order for
ties. Tracks have type `1`; photon and neutral-hadron towers have type `0`.
Tower charge, PID, impact parameters and errors are `-999`, even in occupied
slots. Track values are copied from Delphes; their interpretation depends on
the detector card and configured reconstruction.

| Branch | ROOT type | Units | Description |
|---|---|---|---|
| `num_constituents` | `Int_t` | — | Total number of jet constituents before truncation. |
| `num_constituents_stored` | `Int_t` | — | Number stored: `min(num_constituents, Nc)`. |
| `num_constituents_overflow` | `Int_t` | — | Number omitted because of the capacity. |
| `constituent_pt` | `Float_t[Nc]` | GeV | Track pT or tower transverse energy `ET`. |
| `constituent_eta` | `Float_t[Nc]` | — | Track or tower pseudorapidity. |
| `constituent_phi` | `Float_t[Nc]` | rad | Track or tower azimuth. |
| `constituent_E` | `Float_t[Nc]` | GeV | Track four-vector energy or tower energy. |
| `constituent_type` | `Int_t[Nc]` | — | `1` for a track, `0` for a neutral tower. |
| `constituent_charge` | `Int_t[Nc]` | — | Delphes track charge in units of the elementary charge. |
| `constituent_pid` | `Int_t[Nc]` | — | Delphes track HEP/PDG particle ID. |
| `constituent_d0` | `Float_t[Nc]` | mm | Delphes track transverse impact parameter `D0`. |
| `constituent_dz` | `Float_t[Nc]` | mm | Delphes track longitudinal impact parameter `DZ`. |
| `constituent_d0_error` | `Float_t[Nc]` | mm | Delphes track `ErrorD0`. |
| `constituent_dz_error` | `Float_t[Nc]` | mm | Delphes track `ErrorDZ`. |

### Non-overlapping small-R jets

All remaining Delphes small-R jets are sorted by descending pT; ties retain
their original collection order. They are stored regardless of whether they
individually pass the leptonic-side b-jet requirements.

| Branch | ROOT type | Units | Description |
|---|---|---|---|
| `num_small_r_jets` | `Int_t` | — | Total number of non-overlapping small-R jets before truncation. |
| `num_small_r_jets_stored` | `Int_t` | — | Number stored: `min(num_small_r_jets, Nj)`. |
| `num_small_r_jets_overflow` | `Int_t` | — | Number omitted because of the capacity. |
| `small_r_jet_pt` | `Float_t[Nj]` | GeV | Small-R jet transverse momentum. |
| `small_r_jet_eta` | `Float_t[Nj]` | — | Small-R jet pseudorapidity. |
| `small_r_jet_phi` | `Float_t[Nj]` | rad | Small-R jet azimuth. |
| `small_r_jet_mass` | `Float_t[Nj]` | GeV | Small-R jet mass. |
| `small_r_jet_btag` | `Int_t[Nj]` | — | `1` if Delphes `BTag` bit 0 is set; otherwise `0`. |

For both collections, slots at and beyond the corresponding `_stored` count
are padded with `-999`, including integer arrays. The full count and overflow
are retained when the capacity is exceeded. The b-jet requirement is checked
before small-R truncation: an event can pass even if its qualifying b jet is
not among the stored slots.

### MadGraph truth labels

These values are copied from the matching MadGraph `labels` tree, not derived
from Delphes objects. The current process identifies one final-state bottom
quark, one down-type quark or antiquark (d or s), and one up-type quark or
antiquark by PID and status. The hadronic W is the sum of the two light-quark
four-vectors; the hadronic top is the W plus the bottom. The leptonic top is
the sum of the antibottom, charged lepton and neutrino. See
`madgraph/include/TLorentz_Labels.py` for the calculation.

Four-vectors are actual, unnormalized `(px, py, pz, E)` values in GeV. `lab`
is the input frame, `ttCM` is obtained by boosting to the reconstructed ttbar
center-of-mass frame, and `topRest` is obtained by a further boost from that
frame to the hadronic-top rest frame. Thus the top four-vector in `topRest`
has approximately zero spatial momentum, up to numerical precision. The output
names have the form `<particle>_<component>_<frame>`, with `E` uppercase.

For `costheta_down` and `costheta_bottom`, the reference axis is the unit
direction of the hadronic top in the ttbar CM frame. It is dotted with the
unit direction of the respective daughter after the sequential boost to the
top rest frame. This preserves the existing label convention.

| Branch | ROOT type | Units | Description |
|---|---|---|---|
| `costheta_down` | `Float_t` | — | Cosine for the down-type daughter and the top direction. |
| `costheta_bottom` | `Float_t` | — | Cosine for the bottom daughter and the top direction. |
| `top_{px,py,pz,E}_{lab,ttCM,topRest}` | `Float_t` | GeV | Hadronic top four-vector. |
| `bottom_{px,py,pz,E}_{lab,ttCM,topRest}` | `Float_t` | GeV | Hadronic-side bottom four-vector. |
| `down_{px,py,pz,E}_{lab,ttCM,topRest}` | `Float_t` | GeV | Down-type daughter four-vector. |
| `W_{px,py,pz,E}_{lab,ttCM,topRest}` | `Float_t` | GeV | Hadronic W four-vector. |

Each brace group denotes every combination of one component and one frame:
for example, `top_px_lab` and `W_E_topRest`. There are 12 branches per
particle, plus the two angular-cosine branches.

## Tree metadata

The tree's ROOT user information contains `MAX_LARGE_R_CONSTITUENTS`,
`MAX_SMALL_R_JETS`, `padding_value` and `JET_OVERLAP_DR`. These are metadata,
not branches. Other settings and the radii read from the rendered Delphes card
are recorded in the per-job logs, not in the tree. Only merge successful
outputs with compatible schemas and identical settings.

## Output alignment check

After running the analysis, check the selected jet direction against the
joined MadGraph hadronic-top lab four-vector:

```bash
python3 tests/test_analysis_hadronic_top_largeR_alignment.py \
    --analysis-dir 'analysis/WS_<DELPHES_TAG>'
```

Replace `<DELPHES_TAG>` with the tag used for the output.

This requires `numpy`, `uproot` and `pytest`. It reports the fraction of
selected events with `DeltaR(top, selected large-R jet) <= R`, along with the
unmatched count and distance statistics. By default it requires at least 90%
to match, like the Delphes alignment check. This is a sample-level validation,
not an analysis selection cut or a guarantee that every event matches. The
radius is read from each selected file's analysis log (`R_LARGE`); all logs
must agree. Use `--jet-radius` if logs are unavailable, and
`--min-match-fraction` to change the tolerance. These are validation settings,
not changes to jet reconstruction or event selection.

Alternatively, run:

```bash
python3 -m pytest -q -s tests/test_analysis_hadronic_top_largeR_alignment.py
```

This uses the workspace in `debug.config`. `TP_CONFIG` selects another config,
`TP_ANALYSIS_DIR` overrides the workspace, and `TP_LARGE_R_RADIUS` and
`TP_MIN_MATCH_FRACTION` override the validation settings. The pytest module
also exercises the calculation and its failure cases with synthetic ROOT files,
without regenerating events.
