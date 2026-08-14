# Experiment 22 Amendment A1 — Phase D complete-case index reset

**Trigger:** the once-only formal runner completed Phase M and Phase P, then
crashed in Phase D before writing any Phase D result. The 207-row complete-case
DataFrame retained labels from the 213-row parent table. Gate21's frozen
`within_seed_permutation` implementation correctly assumes a dense positional
index and therefore raised `IndexError: index 207 is out of bounds for axis 0
with size 207`.

**Permitted correction:** reset the complete-case DataFrame index to
`0..n-1` immediately before invoking the unchanged frozen Phase D code.

**Unchanged scientific objects:** all 213 query identities; the 207 evaluable
identities; all predictions, truths and weights; the 110/103 full accepted set;
threshold `0.0923227147328771`; bootstrap/permutation seeds and replicate
counts; all five criteria; geometry control; Phase M/P outputs.

**Reveal status:** no Phase D CSV, NPZ or verdict existed when the defect was
diagnosed. Phase M and Phase P are not rerun. This is an engineering resume
amendment, not a scientific redefinition.

The wrapper monkey-patches only the DataFrame index passed to the original
frozen helper and then executes the original frozen `main()`.
