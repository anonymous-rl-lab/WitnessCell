# Normalization policy

Code 31.0.0 corresponds to paper v41, language revision R1. One root entry and one canonical current manuscript tree replace multiple delivery-layer entry points. Scientific functions and source-specific filenames are retained where fixed records depend on them.

All 180 files of the v41 companion, including its manifest, are represented in provenance/V41_INTEGRATION_MAP.csv. The original source bytes are recoverable exactly. The current physical manifest describes shipped files; ALIASES.json describes byte-identical logical paths; PACKED_ASSETS.json describes lossless gzip payloads; HAN_DERIVED_CACHES.json describes selected weight-derived caches.

Preparation restores compressed payloads first, aliases second, and required Han caches last. These are deterministic input-layout operations, not training. Restored files are excluded from Git so a later commit does not undo compaction. Runtime outputs belong outside the repository.

Older editorial snapshots and obsolete delivery wrappers were removed with SHA256 records. Scientific protocols, negative/positive outcomes and historical comparators required to interpret current claims remain. The root verification directory describes this release; previous-run records are provenance.
