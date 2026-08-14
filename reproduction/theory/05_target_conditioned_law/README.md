# Held-out target-conditioned geometry diagnosis

This is the first CPU-only gate after the original topology demonstration.
It fixes the prediction-evaluation weakness of the earlier probe:

- every evaluated target pair is absent from training;
- every target has both endpoints represented in training (`seen2`);
- train-edge count and global nullity are fixed;
- witnessed and unwitnessed target directions coexist in the same design;
- outcomes are never used to compute the structural risk score.

Run a smoke test:

```bash
python run_target_geometry_diagnosis.py --seeds 5 --out /tmp/witness_smoke
```

Run the frozen CPU test:

```bash
python run_target_geometry_diagnosis.py --seeds 100
```

The preregistered gate requires: protocol audit pass; geometry failure AUROC
at least 0.90; at least 0.10 AUROC above the pair-type-only baseline; and at
least a three-fold NMSE separation between unsupported directions and all
other targets.

Passing this test validates the algebraic H1 in a controlled system.  It does
not yet validate target-conditioned acquisition on Norman or improve an AIVC.

## Outcome-blind acquisition gate

After the diagnosis passes, compare equal-budget acquisition policies:

```bash
python run_acquisition_trial.py --seeds 5 --random-reps 5 --out /tmp/acquisition_smoke
python run_acquisition_trial.py --seeds 50 --random-reps 20
```

The candidate and target outcomes are hidden during selection.  Global
D-optimal and target-conditioned V-optimal obtain the same global-nullity
reduction; V-optimal alone may use the predeclared target identities.  The
primary budget equals the number of target-bearing components.

## Geometry–adequacy phase diagram

```bash
python run_adequacy_geometry_phase.py --seeds 3 --out /tmp/phase_smoke
python run_adequacy_geometry_phase.py --seeds 20
```

This adds heterogeneous pair-specific mismatch while keeping all targets
strictly held out and seen2. A training-only OOF adequacy gate routes risk
assessment to geometry in factorized regimes and to endpoint representation
bias when the factorization breaks.
