# Witness geometry synthetic truth test

This experiment isolates the exact graph-theoretic boundary for the anchor
model `R_ij = a_i + a_j`.  A dense bipartite intervention graph leaves one
unidentifiable contrast even though every node has many observed partners.  A
single odd-cycle witness removes that null direction.

```bash
python simulate_topology_law.py --seeds 20 --out results
```

The experiment is a planted identifiability test, not biological evidence.
