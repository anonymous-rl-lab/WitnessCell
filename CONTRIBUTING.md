# Contributing

Keep scientific results traceable to permitted inputs, fixed splits and the named task. Submit code or documentation changes separately from changes to a frozen analysis. Do not overwrite reference predictions or expected results to make a check pass. Add a new analysis identity for changed estimation rules.

Run `python -B -m unittest discover -s tests -v` and the affected reproduction task. Record the environment and numerical differences. Regenerate release manifests only when intentionally building a new version, with changes reviewed against the previous manifest. Local output should remain outside the repository.

Public contribution and distribution arrangements depend on the licence scope described in LICENSE.md.
