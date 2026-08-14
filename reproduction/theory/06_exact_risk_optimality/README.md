# Exact risk and optimality probes

This experiment replaces fixed score averaging and hard routing with the exact
best-linear-unbiased predictive risk derived in `THEORY.md`.

Order of execution:

```bash
python verify_exact_risk.py --replicates 20000 --random-competitors 100
python enumerate_optimal_design.py --instances 3 --mc-replicates 5000 --out /tmp/design_smoke
python enumerate_optimal_design.py --instances 30 --mc-replicates 10000
```

The first script must pass before any design comparison is interpreted. The
enumeration experiment uses known covariance and therefore tests oracle
optimality. Estimating that covariance from real perturbation data is outside
this gate.
