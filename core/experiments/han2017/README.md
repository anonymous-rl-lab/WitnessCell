# Han2017 compact scientific snapshot

The paper-relevant protocol, all configurations, original primary and later full-screen extension are preserved. Read `../../docs/HAN2017_PAPER_SUMMARY_CN.md` .

`snapshot/` contains unchanged scientific files. Its historical manifest describes the original archive, including deterministic files now described in `RECIPES.json`; use the root release manifest and `COMPACTION_LEDGER.json` for this compact state. Run `materialize.py --out <fresh external path> --mode exact` to restore the original array layout. The full validated entry is root `tools/replay_han.py`; direct old archive paths are historical.
