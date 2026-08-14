# Third-party notices

The WitnessCell Python source distribution does not vendor third-party source
code, datasets, GO annotations, pretrained models, or manuscript files.
Runtime and development dependencies are installed separately and retain their
respective licenses and notices.

The package is distributed under Apache-2.0. Dependency declarations are in
`pyproject.toml`; downstream distributors are responsible for reviewing the
licenses of the exact dependency versions they redistribute.

The GitHub reproducibility layer additionally contains derived condition
moments, GO assets, frozen benchmark tables and independent interoperability
implementations. Those files are not included in the PyPI wheel or sdist.

Experiment 22 is source-locked to the MIT-licensed
`shiftbioscience/Perturbation-Models-Outperform-Baselines` repository at commit
`b2ff09f5f83a6e7a31a13d25d046f504f7325178`. Comparator semantics are
source-locked to the GPL-3.0 `bm2-lab/scPerturBench` repository at commit
`6e24e7a9827e55d4567d2139427be9af0d1e7a6c`; its source is not vendored here.
The detailed provenance and redistribution boundary is preserved in
`reproduction/gate22_metrics/THIRD_PARTY_NOTICE.md`.

Public single-cell matrices are downloaded separately from their official
release and retain the original data terms. Users redistributing upstream data
or source must retain the corresponding licenses and notices.
