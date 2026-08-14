# Gate 19 — WitnessCell v14 incremental amplitude gate

Gate 19 is the formal v14 integration layer on top of the unchanged v13 release.
It preserves the v13 endpoint background, factorized backbone, interaction
witness kernel, official GEARS splits, validation selection and scoring protocol.

For an unseen endpoint, v14 adds a one-coordinate amplitude-residual correction
estimated only from training singles. Honest outer leave-one-endpoint-out fits
activate the correction only when both one-sided 95% lower bounds are positive:

1. all-gene MSE improvement over a v13 refit;
2. top-100 PCC improvement over the same v13 refit.

The gate activates in 8/12 official splits. Wessels seeds 1–3 and Schmidt seed 3
return v13 predictions exactly. Across 654 condition×seed units, v14 improves
all-gene MSE by 0.108%, top-100 MSE by 0.305%, and top-100 PCC by 0.00130
relative to v13.

The release contains 12 target-free deployment archives. The directional
top-100/top-5000 comparison is rerun and remains rank 1. The official
distributional six-metric tables under Gate 18 are inherited historical v13
evidence; they are not relabeled as newly scored v14 results.

See `ALGORITHM_CONTRACT.md`, `FORMAL_REPORT_CN.md`, `RUNBOOK.md`, and
`RELEASE_AUDIT.json`.
