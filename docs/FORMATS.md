# File formats

## Condition moments (`.npz`)

`ConditionMoments.from_npz` requires a NumPy archive loadable with
`allow_pickle=False` and these exact logical arrays:

| Key | Shape | Type | Meaning |
| --- | --- | --- | --- |
| `genes` | `(G,)` | Unicode | Unique gene order |
| `conditions` | `(C,)` | Unicode | Unique condition labels |
| `means` | `(C, G)` | numeric | Condition expression means |
| `variances` | `(C, G)` | numeric | Population variances (`ddof=0`) |
| `counts` | `(C,)` | integer | Positive cell counts |

Use `ConditionMoments.to_npz` to produce a conforming archive.

## Split JSON

```json
{
  "train_conditions": ["control", "A", "B", "C", "A+B"],
  "validation_conditions": ["A+C"]
}
```

The arrays must be disjoint. Final prediction targets must not be listed.

## Gene-to-GO JSON

```json
{
  "A": ["GO:0000001", "GO:0000002"],
  "B": ["GO:0000002"]
}
```

Keys and terms are normalized to strings. No executable pickle input is used by
the public CLI.

## Model bundle (`.wcell`)

A model is a ZIP with exactly three files:

- `manifest.json`: format/version plus SHA-256 and size for both payloads;
- `metadata.json`: JSON configuration, gates, diagnostics and named state; and
- `arrays.npz`: numeric arrays loadable with `allow_pickle=False`.

Entries are deterministic and have fixed timestamps. Saving is atomic within a
filesystem. Loading rejects unexpected/duplicate entries, duplicate JSON keys,
unsupported format versions, oversize or over-expanded nested archives,
hash/size mismatches, non-float64 or non-finite arrays, invalid fitted ranges,
configuration/state mismatches, leakage overlaps and inconsistent shapes.
Integrity is not a digital signature; authenticate distribution files with the
release checksums or publisher attestations.

## Predictions (`.npz`)

`PredictionBatch.to_npz` writes Unicode/numeric arrays only:
`conditions`, `genes`, `means`, `effects`, `factorized_means`, and
`factorized_effects`, plus optional `risks` and `decisions`. It never stores
target truth.
