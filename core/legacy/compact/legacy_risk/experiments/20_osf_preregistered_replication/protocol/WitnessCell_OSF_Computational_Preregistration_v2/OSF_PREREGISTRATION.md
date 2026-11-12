# OSF preregistration

## 1. Study information

### Title

Preregistered cross-environment computational reproduction of WitnessCell v14

### Study type

This is a confirmatory computational reproducibility study. It regenerates
frozen WitnessCell v14 outputs from packaged training-side condition-moment
caches and audits the released aggregation artifacts. It involves no new cell
culture, perturbation, sequencing, guide construction, biosafety assessment,
human participants, animal work or biological sample collection.

### Scientific question

Can an independent Python 3.12 execution reproduce the frozen WitnessCell v14
deployment outputs, evidence-gate decisions and paired mean-response
comparisons, while also verifying the integrity of the explicitly inherited
Gate18 scPerturBench aggregation, without modifying the algorithm, split, cache
or scoring contracts?

### Relation to earlier work

The earlier OSF registration at `https://osf.io/a2dfr` specified an unperformed
wet-lab extension and has been withdrawn. No prospective wet-lab outcome was
collected. Before this v2 registration, one exploratory AutoDL execution was
observed to reproduce the frozen computational outputs. That prior result is
fully disclosed and is not counted as the registered confirmatory run. The
present registration freezes the protocol for one new execution performed only
after OSF accepts this v2 registration.

The purpose is confirmation of computational reproducibility, not discovery of
an unknown biological effect. Exact expected values are therefore specified in
advance.

## 2. Frozen inputs and exclusions

### Included input

The only source archive is `grn_repo_v14.zip`, SHA-256
`955d446a0cfe65a6f4f191d74c7a7d16a6083488b9beb40a5e155ab94d528bd7`.
The archive contains the frozen v14 code, v13 comparator archives, official
split manifests, training-side condition-moment caches, target-free v14
deployment archives and packaged benchmark aggregation records.

### Excluded inputs

The registered execution will not use:

- new biological measurements;
- the withdrawn 135-pair panel or its former laboratory documents;
- raw prospective target outcomes;
- any post-registration algorithm, split, cache or threshold change;
- any result from the earlier AutoDL execution as an execution input;
- raw h5ad rescoring of the six distributional benchmark metrics.

The regenerated `predictions.npz` files contain historical held-out benchmark
truth solely for the frozen paired audit. Only `deploy_predictions.npz` is a
target-free deployable artifact.

## 3. Frozen execution design

### Environment

The confirmatory run uses Python 3.12 with:

- NumPy 1.26.4
- pandas 2.2.3
- SciPy 1.14.1
- scikit-learn 1.5.2

The transitive numerical dependencies are also pinned in `requirements.txt`
and verified by the execution entry point.

The run is CPU-only. No model is retrained on new observations. Four datasets
(`Norman`, `Replogle_exp6`, `Schmidt`, `Wessels`) and seeds 1, 2 and 3 produce
12 deterministic split executions.

### Frozen command

The sole confirmatory entry point is:

```bash
bash WitnessCell_OSF_Computational_Preregistration_v2/run_preregistered_replication.sh /root/autodl-tmp/grn_repo_v14.zip
```

The script verifies the source hash and environment, creates a new output
directory, regenerates all 12 splits, runs the paired v14-versus-v13 audit and
applies the fail-closed registered audit. Manual replacement of an output,
rerunning only failed splits after inspecting their numerical results or
editing the verdict is forbidden.

## 4. Confirmatory criteria

### C1 -- execution completeness and schema

Exactly 12 regenerated target-free deployment archives must exist, one for
each dataset--seed split. Each archive must contain exactly the same 13 named
arrays, shapes and dtypes as its frozen v14 counterpart, and contain neither
`truth` nor `truth_variance`.

### C2 -- numerical identity

All 156 regenerated deployment arrays (13 arrays x 12 splits) must be
elementwise identical to the corresponding frozen v14 arrays. The global
maximum absolute difference over numeric arrays must equal 0. String and
categorical arrays must also be exactly equal.

### C3 -- evidence-gate and fallback identity

The regenerated manifests must recover eight correction-active splits and four
inactive exact fallbacks. The inactive set must be exactly Schmidt seed 3 and
Wessels seeds 1--3. Every inactive prediction must be exactly equal to the
frozen v13 prediction.

### C4 -- paired mean-response reproduction

Across exactly 654 condition--seed units, the regenerated v14-versus-v13 audit
must reproduce, with absolute tolerance `1e-12`:

- relative all-gene MSE improvement: `0.001079217249737252`;
- relative top-100 MSE improvement: `0.0030453806218846413`;
- top-100 Pearson-correlation change: `0.0013006984799760163`.

### C5 -- inherited Gate18 benchmark-aggregation integrity

Gate19/v14 did not recompute the six distributional metrics from raw h5ad. Its
release explicitly inherited the Gate18 benchmark evidence. As a secondary
artifact-integrity criterion, that frozen Gate18 official top-100 aggregation
must contain 654 operation--seed units and reproduce:

- WitnessCell mean operation rank: `3.871559633027523`;
- official position: 1;
- both-panel robustness position: 1;
- dataset-equal top-100 position: 1;
- scouter comparison rank difference: `-1.5412844036697249`;
- scouter dataset--perturbation cluster-bootstrap 95% interval:
  `[-1.9540636702850922, -1.1267079136021945]`.

C5 audits inherited frozen benchmark tables. It is not a claim that the
confirmatory run recomputed Wasserstein, energy distance, KL or common-DEG
scores from raw h5ad matrices, and it must not be described as a new v14
six-metric result.

## 5. Success, failure and stopping rules

The registered confirmatory result is `PASS` only if C1, C2, C3, C4 and C5 all
pass in the same uninterrupted execution and the source hash matches. Any
missing split, package mismatch, environment mismatch, non-finite output,
schema change, nonzero array difference or criterion failure yields `FAIL`.

An infrastructure interruption before a verdict is `NOT_EVALUABLE`. One clean
restart from the beginning is permitted only when no criterion-level numerical
result was produced. Partial output cannot be substituted. A failed numerical
criterion cannot be repaired by changing code or thresholds under this
registration.

## 6. Reporting plan

The result archive, environment record, full log, generated verdict and SHA-256
digest will be uploaded unchanged. All criteria will be reported, including
failures. The prior exploratory execution and the withdrawn earlier
registration will remain disclosed.

If all criteria pass, the manuscript may state that the frozen WitnessCell v14
computational object was reproduced in a preregistered confirmatory
cross-environment execution. It may not describe this study as prospective
biological validation, a new perturbation experiment or confirmation of the
withdrawn H1--H5 biological hypotheses.
