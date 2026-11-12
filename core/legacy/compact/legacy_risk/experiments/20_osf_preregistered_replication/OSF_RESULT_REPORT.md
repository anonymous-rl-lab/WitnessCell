# Result report: WitnessCell v14 confirmatory computational replication

## Registration reference

- OSF registration: `xrsmd`
- URL: https://osf.io/xrsmd
- Protocol version: `2.0.0`
- Study type: preregistered confirmatory computational replication
- Prior successful execution disclosed in the registered protocol: yes

## Frozen inputs

- WitnessCell source version: `14.0.0`
- Source archive SHA-256: `955d446a0cfe65a6f4f191d74c7a7d16a6083488b9beb40a5e155ab94d528bd7`
- Preregistration package SHA-256: `8d2bcc1c06a42fee039a515ee12588d516eb04c9979697181eb35ebd4c4a8c63`

## Execution

- Platform: AutoDL, CPU-only registered workflow
- Start: `2026-08-10T10:56:15Z`
- Completion: `2026-08-10T11:13:03Z`
- Duration: 1,008 seconds (16 minutes 48 seconds)
- Python: `3.12.13`
- NumPy: `1.26.4`
- pandas: `2.2.3`
- SciPy: `1.14.1`
- scikit-learn: `1.5.2`
- joblib: `1.5.3`
- threadpoolctl: `3.6.0`

## Prespecified criteria

| Criterion | Frozen requirement | Result |
|---|---|---|
| C1 | 12 deployment archives with 13 matching target-free arrays per archive | **PASS**: 12 archives and 156 array comparisons |
| C2 | 156/156 arrays elementwise exact; global maximum absolute difference 0 | **PASS**: 156/156 exact; maximum difference `0.0` |
| C3 | Eight active splits and four exact fallbacks | **PASS**: eight active; exact fallback for Schmidt seed 3 and Wessels seeds 1–3 |
| C4 | Exact recovery of the frozen 654-unit paired mean-response summary | **PASS**: all-gene MSE gain `0.001079217249737252`; top-100 MSE gain `0.0030453806218846413`; top-100 Pearson delta `0.0013006984799760163` |
| C5 | Integrity of inherited Gate18 official aggregation | **PASS**: mean operation rank `3.871559633027523`; position 1; both-panel position 1; dataset-equal top-100 position 1 |

For C5, the frozen WitnessCell-minus-scouter rank difference was `-1.5412844036697249`, with dataset–perturbation cluster-bootstrap 95% interval `[-1.9540636702850922, -1.1267079136021945]`.

## Overall decision

`PASS_PREREGISTERED_COMPUTATIONAL_REPLICATION`

The post-registration run exactly reproduced the released WitnessCell v14 deployment object across a second software environment. It also reproduced the prespecified v14-versus-v13 paired summary and verified the integrity of the explicitly inherited Gate18 aggregation.

## Returned result archive

- Filename: `WitnessCell_preregistered_replication_20260810T105615Z.tar.gz`
- Size: `39,357,228` bytes
- SHA-256: `511860f645158a9af536a3944e03d4cc440674e23b55d826bf615e02d96bea0c`

## Interpretation boundary

This result establishes exact computational reproducibility under the registered protocol. It does not constitute a new biological measurement, a wet-lab validation, a prospective biological outcome reveal or a new Gate19/v14 raw-h5ad six-metric evaluation. C5 is an integrity audit of packaged Gate18 evidence.

