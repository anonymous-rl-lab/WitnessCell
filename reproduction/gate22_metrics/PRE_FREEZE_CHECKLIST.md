# Experiment 22 v2 pre-freeze checklist

Scientific freeze and formal execution are blocked until every required row is `PASS`. `NOT_APPLICABLE` is allowed only where stated. A missing asset produces `NOT_ADJUDICATED`, never an implicit pass.

| ID | Requirement | Evidence required | Current state |
|---|---|---|---|
| S1 | normative repository resolves to the locked commit/tree | `source_audit.json` SHA `7f6ee536…614670` | PASS |
| S2 | all locked source-file SHA-256 values match | `source_audit.json` SHA `7f6ee536…614670` | PASS |
| S3 | upstream environment is reproducible under Python ≥3.12 | `environment_report.json` SHA `0290da80…b9e879`; smoke/source tests | PASS |
| S4 | local metric wrapper is source-parity tested | 16 tests, `unit_tests.log` SHA `34ff20f1…a7580` | PASS |
| S5 | zero-weight, duplicate-gene, NaN and reordered-gene cases fail or resolve exactly as specified | negative tests in `unit_tests.log` | PASS |
| A1 | four public formal cell-level datasets are present with accession/checksum | `cell_asset_ledger.json` SHA `352c49d3…e771d` | PASS |
| A2 | official 12 split identities and 654 unit identities match Gate 19 | `frozen_asset_audit.json` SHA `7c47f9f1…ace3c8e` | PASS |
| A3 | technical-duplicate cell-ID splits are generated without prediction-file access | `CONTRACT_SEAL.json` SHA `a963aa5a…30b9c`; four firewall-clean runners | PASS |
| A4 | truth-side DEG weights reproduce locked cell-level code | `weight_contract_audit.json` SHA `02186a5f…357b` | PASS |
| A5 | v14 and v13 prediction archives pass their existing manifests | `frozen_asset_audit.json` | PASS |
| C1 | trainMean and baseControl are exactly regenerated | comparator source audit + branch unit test; audit SHA `c9a98728…194ad` | PASS |
| C2 | linearModel availability is adjudicated before scoring; absent raw predictions/runtime force `PRED_LINEAR=NOT_ADJUDICATED` | `comparator_audit.json`; status fixed before scoring | PASS |
| C3 | remaining comparator raw means are inventoried individually | `comparator_inventory.csv` SHA `c41a881c…377de` | PASS |
| C4 | comparator scope label is fixed before prediction scoring | `COMPARATOR_SCOPE.json` SHA `46b0a193…3a2f8` | PASS |
| G1 | Gate 19 active/inactive labels match 8/4 frozen split ledger | exact identity audit | PASS |
| G2 | four inactive v14 arrays are exactly equal to v13 arrays | exact array audit | PASS |
| G3 | shadow-gate pseudo-target weights use training LOO only | `shadow_source_flow.json` SHA `fdeff1a3…6f06d` | PASS |
| Q1 | Gate 21 threshold is bitwise `0.0923227147328771` | frozen asset + gene-contract audits | PASS |
| Q2 | Gate 21 accepted query IDs exactly match the original result | accepted-ID SHA `3f6c0c7c…55d75` | PASS |
| Q3 | Gate 21 gene-level predictions and truths are regenerated | gene contract SHA `d58cec87…4c41`; report SHA `cd6689e8…e8b6` | PASS |
| Q4 | Gate 21 canonical weight availability is adjudicated for all 213 identities; missing rows force full-panel `NOT_ADJUDICATED` and remain only in a labelled ≥90% complete-case sensitivity | 207/213 rows; missing pair fixed as `ELMSAN1+MAP2K3`; no imputation | PASS |
| F1 | Phase M environment cannot read candidate predictions | intentional denial + candidate-blind smoke SHA `ca0e862b…22c7e` | PASS |
| F2 | Phase P cannot alter Phase M contract | `phase_boundary_audit.json` SHA `8e373ec4…14f` | PASS |
| F3 | all formal output paths are empty before freeze | phase-boundary audit: `formal_e22` absent | PASS |
| F4 | protocol, code, assets and environment are hashed | `FROZEN_MANIFEST.sha256`; external hash in `FREEZE_RECEIPT.json` | PASS |
| D1 | superseded v1 draft is archived and excluded from runners | `SUPERSESSION_NOTICE.md` SHA `592fce74…19f2eee` + runner audit | PASS |

## Required source-parity tests

- perfect prediction: WMSE = 0; weighted delta R² = 1; NIR = 1 when every matched truth is uniquely nearest;
- source fixed arrays: local values equal locked functions at the stored dtype/tolerance;
- permutation: NIR changes according to the locked strict comparison, including ties;
- DRF: lower- and higher-is-better branches equal the locked implementation, including epsilon and clipping;
- zero-sum weights: formal unit marked non-evaluable; no uniform fallback;
- gene order: reordering predictions without reordering names fails closed;
- duplicate genes: maximum locked weight retained before final alignment;
- phase firewall: Phase M raises on every prediction path.

## Freeze action

After all rows pass:

1. replace `Current state` with the evidence status and exact artifact hash;
2. create the implementation and asset manifest;
3. change protocol status to `FROZEN_BEFORE_EXPERIMENT22_FORMAL_EXECUTION`;
4. generate `FROZEN_MANIFEST.sha256` once;
5. do not edit a frozen file—create a versioned amendment if a defect is found.
