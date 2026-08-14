# Official repository integration

1. Apply `calPerformance_genetic.patch` at the root of scPerturBench.
2. Keep the bundled `data/gears_assets` directory beside the adapter, or set
   `WITNESSCELL_GEARS_ASSETS` to an existing verified copy. These files
   reproduce the GO-support filter applied before the official GEARS split;
   they are read-only and do not install GEARS.
3. Copy `myWitnessCell.py` into
   `Perturbation_generalization/Genetic/myWitnessCell.py`.
4. Export files to the standard location:
   `DataSet2/<dataset>/hvg5000/WitnessCell/savedModels<seed>/result.h5ad`.

The local audit used scPerturBench commit
`6e24e7a9827e55d4567d2139427be9af0d1e7a6c`.
