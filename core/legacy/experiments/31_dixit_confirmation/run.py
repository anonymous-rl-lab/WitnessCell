#!/usr/bin/env python3
"""One entry: primary RNA summary audit and distinct G2M development refit."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if out == ROOT or ROOT in out.parents:
        parser.error('Use an output directory outside the immutable experiment module.')
    if out.exists() and any(out.iterdir()):
        parser.error('Use a new or empty output directory.')
    out.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, str(ROOT / 'replay_primary_summary.py'),
                    '--out', str(out / 'PRIMARY_RNA_SUMMARY_AUDIT.json')], check=True)
    subprocess.run([sys.executable, str(ROOT / 'REPRODUCE.py'),
                    '--output', str(out / 'g2m_development')], check=True)
    primary = json.loads((out / 'PRIMARY_RNA_SUMMARY_AUDIT.json').read_text())
    development = json.loads((out / 'g2m_development/REPRODUCTION_COMPLETE.json').read_text())
    result = {'status': 'PASS', 'primary_rna': primary,
              'g2m_development_named_arrays_exact': development['exact_named_arrays'],
              'g2m_is_independent_confirmation': False,
              'primary_rna_predictions_refitted': False,
              'g2m_development_refitted_from_condition_means': True}
    (out / 'INTEGRATED_REPLAY.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
