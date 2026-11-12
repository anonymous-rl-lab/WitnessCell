# Gate 20: OSF-preregistered cross-environment replication

This gate adds the completed result of the confirmatory computational-replication protocol registered at [OSF xrsmd](https://osf.io/xrsmd) to the WitnessCell release tree.

## Release role

Gate 20 does **not** change the WitnessCell v14 numerical algorithm, frozen predictions or benchmark scores. It adds a second-environment reproducibility layer:

1. the exact v2 preregistration package and its expanded files;
2. the original AutoDL result archive, unchanged;
3. extracted C1–C5 evidence and target-level summaries;
4. an independent audit that compares all regenerated deployment arrays directly with Gate 19;
5. explicit claim boundaries for the OSF result.

## Frozen identities

| Asset | SHA-256 |
|---|---|
| v14 source archive used by the registered run | `955d446a0cfe65a6f4f191d74c7a7d16a6083488b9beb40a5e155ab94d528bd7` |
| v2 preregistration package | `8d2bcc1c06a42fee039a515ee12588d516eb04c9979697181eb35ebd4c4a8c63` |
| returned AutoDL result archive | `511860f645158a9af536a3944e03d4cc440674e23b55d826bf615e02d96bea0c` |
| OSF result-registration package | `0fa75b22fd698a2e72350ba3ffd8eb680da6c94ada69cb19f12ca98c7927bd26` |

## Confirmatory result

`PASS_PREREGISTERED_COMPUTATIONAL_REPLICATION`

- 12/12 target-free deployment archives regenerated;
- 156/156 arrays elementwise exact;
- global maximum absolute difference `0.0`;
- eight active correction splits and four exact fallbacks;
- 654 paired condition–seed units reproduced;
- C1–C5 passed in one post-registration execution.

Run the independent release audit from the repository root:

```bash
python experiments/20_osf_preregistered_replication/audit.py
```

## Claim boundary

This is a preregistered confirmatory computational replication following a disclosed earlier successful execution. It is not a wet-lab experiment, new biological data collection, prospective biological outcome reveal or new Gate19/v14 raw-h5ad six-metric rescoring. Criterion C5 checks the integrity of the inherited Gate18 aggregation.

