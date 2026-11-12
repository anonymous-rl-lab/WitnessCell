# Canonical manuscript: v41 R1

`manuscript_sources/` contains the delivered manuscript, SI, captions and vectors. The manuscript records its original two-archive Code availability; both components are integrated in this code release, so separate ZIP files are no longer runtime dependencies. The root README is the current operational guide.

The older figure generator in `code/` is a dependency for unchanged panels; the root figure task applies the current reference-sensitivity generator and checks all 12 panels against `reference_figures/`. Table generators and the root coverage map distinguish calculated values from descriptive source/analysis-status rows.

Use `python -B reproduce.py --out ../paper-pdf --tasks paper_pdf` from the repository root to render PDFs. Historical filenames under core preserve fixed source relationships and do not identify the current manuscript.
