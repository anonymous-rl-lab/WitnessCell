# Third-party source notice

Experiment 22's metric definitions are source-locked to
`shiftbioscience/Perturbation-Models-Outperform-Baselines` at commit
`b2ff09f5f83a6e7a31a13d25d046f504f7325178`. That repository identifies its
original portions as MIT-licensed (Copyright 2025 Shift Bioscience). The local
`metric_core.py` is a small, fail-closed interoperability implementation of the
specified metric arithmetic; source-parity tests execute against the locked
upstream functions.

Comparator semantics are source-locked to `bm2-lab/scPerturBench` at commit
`6e24e7a9827e55d4567d2139427be9af0d1e7a6c`, whose repository license is GPL-3.0.
Experiment 22 does not vendor those source files. `comparator_core.py` provides
an independently readable Python execution of the condition-mean equations
exposed by the released scripts, with immutable file hashes and attribution in
`COMPARATOR_SOURCE_LOCK.json`.

No third-party model predictions are reconstructed from published summary
scores. Users redistributing upstream files must retain their original license
and notices.
