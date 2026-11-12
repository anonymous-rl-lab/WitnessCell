"""Unified repository entry for the compact, read-only SF12 experiment replay."""
import sys
sys.dont_write_bytecode = True
import argparse
import json
from run_replay import replay

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, help="External directory for generated QA reports")
    a = p.parse_args()
    result = replay(a.out)
    print(json.dumps({k: result[k] for k in ["status", "score_max_abs", "SF12_coefficient_max_abs", "case_R_coefficient_max_abs", "case_R_quadratic_max_abs", "contrast_identity_rank", "unique_query_edges"]}, indent=2))
