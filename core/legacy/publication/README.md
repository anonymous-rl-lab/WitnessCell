# Paper-v22 publication bridge

This compact module keeps the four original figure builders and only small derived
source tables / signed contrasts used by the published figures and supplementary
tables. Original builders and copied data are byte-identical to the final v22
paper package. No manuscript PDFs, old versions, training matrices or checkpoints
are stored here.

Run `python run.py --out NEW_EXTERNAL_DIRECTORY` to verify source identities and
rebuild 7 main figures plus 10 Extended Data figures (PDF, PNG and SVG). The
historical Fig2/Fig3 remapping is applied automatically. Only Python, NumPy and
Matplotlib are required. Outputs and copied plotting inputs stay outside this
immutable module.

The builder also independently checks the saved 11×4,995 RNA contrast errors,
four-corner example sums, 1,234 Jurkat repeat rows, benchmark counts and the top-ten
gain contribution. These are saved-output checks. Model fits, preprocessing and
selection are verified by the corresponding experiment modules.

`PUBLISHED_TABLES_v22.md` is the exact typeset reference for all 19 supplementary
tables. It is an editorial comparison target, not a substitute for underlying
reproduction. `../docs/PAPER_COVERAGE.csv` maps every figure/table to its experiment and
reproduction depth. Original-system extension ED9/S11/S12 is reproducible at the
archived split-table aggregation level; this module does not invent its omitted
full predictions or claim raw-model retraining. The runner recomputes all 54
split-method gains and 18 three-seed dataset aggregates from these archived
loss tables, using their original equal-seed geometric aggregation.

The explicit numerical data here may overlap with small presentation exports
inside unchanged compact E25–32 archives. Those experiment archives are preserved
unchanged at the user's request; no additional raw inputs are introduced.
