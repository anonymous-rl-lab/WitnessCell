# Numerical fidelity

1. **File integrity:** SHA256 and size are exact for shipped assets and losslessly restored aliases. Input numeric arrays are not rounded, quantized or converted to a lower precision.
2. **Han cache restoration:** strict mode requires all 120 reconstructed INNER/QUERY NPZ caches to match saved bytes. Compatibility mode records their hashes and leaves numerical acceptance to downstream score/weight checks; it is not a relaxed pass for scientific values.
3. **Scoring:** each scientific task retains its original explicit tolerances. Current reference-sensitivity replay checks nine tables at `rtol=atol=1e-10`; its independent direct-design solves require maximum difference below `1e-9`.
4. **Figures:** data-derived vectors are regenerated; rasterized comparison excludes metadata timestamps. The current figure task checks all 12 figures at the previously defined pixel thresholds (mean absolute difference below 0.05; changed channel fraction below 0.002).
5. **PDFs:** manuscript content is versioned; binary-identical PDFs across TeX/font environments are not promised.

The released validation reports state observed differences. Historical validation records are stored under provenance and are not presented as new runs. A successful score replay does not establish that FASTQ processing or all predictor training was repeated.
