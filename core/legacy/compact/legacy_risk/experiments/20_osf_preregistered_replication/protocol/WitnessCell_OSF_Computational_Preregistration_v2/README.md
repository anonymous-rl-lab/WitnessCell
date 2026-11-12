# WitnessCell computational preregistration v2

This package preregisters a confirmatory, cross-environment computational
reproduction of the frozen WitnessCell v14 release. It contains no wet-lab
study, no new biological sample collection and no prospective biological
outcome claim.

## Registration order

1. Upload this complete package to a **new** OSF registration.
2. Record the new immutable registration URL and timestamp.
3. Do not edit the source archive or this protocol after registration.
4. Only after registration, run the command below in a clean Python 3.12
   environment.
5. Upload the returned result archive and its SHA-256 digest as the registered
   outcome.

The withdrawn earlier registration at `https://osf.io/a2dfr` is not replaced or
silently rewritten. Its scope and withdrawal are disclosed in
`PRIOR_EXECUTION_DISCLOSURE.md`.

## Frozen source

- Source archive: `grn_repo_v14.zip` (the uploaded file may be renamed)
- Required SHA-256:
  `955d446a0cfe65a6f4f191d74c7a7d16a6083488b9beb40a5e155ab94d528bd7`
- WitnessCell version: `14.0.0`

## One-command execution

After extracting this registration package and placing the v14 source archive
on AutoDL, run:

```bash
bash WitnessCell_OSF_Computational_Preregistration_v2/run_preregistered_replication.sh /root/autodl-tmp/grn_repo_v14.zip
```

If the source archive has another path or name, replace only the final path.
The script verifies its SHA-256 before execution and fails closed on mismatch.

## Primary registered claim

If and only if every frozen criterion C1--C5 passes, the result may be described
as:

> A preregistered confirmatory computational replication, conducted after the
> v2 registration and following a disclosed prior successful execution,
> reproduced the frozen WitnessCell v14 prediction archives, evidence-gate
> branches and paired mean-response results, together with the integrity of the
> explicitly inherited Gate18 benchmark aggregation.

This wording does not claim new biological validation, Gate19/v14 raw-h5ad
distributional rescoring or a newly computed v14 six-metric table.
