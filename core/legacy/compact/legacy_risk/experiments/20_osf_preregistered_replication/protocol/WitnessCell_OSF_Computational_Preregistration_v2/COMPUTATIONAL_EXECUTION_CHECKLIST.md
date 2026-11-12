# Computational execution checklist

## Before registration

- [ ] Confirm the earlier registration is visibly marked withdrawn.
- [ ] Upload this complete v2 package to a new OSF registration.
- [ ] Confirm the v2 registration is immutable and record its URL and time.
- [ ] Do not execute the confirmatory command before that timestamp.

## Before execution

- [ ] Place the frozen v14 source archive on AutoDL.
- [ ] Confirm no earlier result directory is being reused.
- [ ] Use Python 3.12 with the versions in `requirements.txt`.
- [ ] Do not edit code, caches, split manifests or expected values.

## Execute once

```bash
bash WitnessCell_OSF_Computational_Preregistration_v2/run_preregistered_replication.sh /root/autodl-tmp/grn_repo_v14.zip
```

## After execution

- [ ] Record whether the status is `PASS`, `FAIL` or `NOT_EVALUABLE`.
- [ ] Do not rerun a numerical failure under the same registration.
- [ ] Upload the complete returned `.tar.gz` result archive.
- [ ] Upload or record its SHA-256 digest.
- [ ] Report C1--C5 individually.
