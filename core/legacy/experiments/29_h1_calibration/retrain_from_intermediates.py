#!/usr/bin/env python3
"""Optional historical full-expression rerun; upstream arrays are not bundled.

Uses the unchanged original trainer after resolving paths explicitly and
verifying the four upstream files. This is not the compact default replay.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', required=True, type=Path,
                   help='Original batch_preprocessor_20261020 directory containing inputs/ and stage2_inputs/.')
    p.add_argument('--output', required=True, type=Path,
                   help='Empty or new output directory; never the reference directory.')
    args = p.parse_args()
    inputs = args.inputs.resolve(); output = args.output.resolve()
    expected = json.loads((HERE / 'reference/SOURCE_INPUT_HASHES.json').read_text())
    for rel, digest in expected.items():
        path = inputs.joinpath(*Path(rel).parts[1:])
        if not path.is_file():
            raise FileNotFoundError(f'Required upstream intermediate is absent: {path}')
        h = hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda: f.read(1 << 20), b''):
                h.update(block)
        if h.hexdigest() != digest:
            raise RuntimeError(f'Upstream source hash mismatch: {path.name}')
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('Output must be empty; archived outputs are not overwritten.')
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(HERE / 'source/origin_free_model_frozen.py', output / 'origin_free_model_frozen.py')
    sys.path.insert(0, str(HERE / 'source'))
    spec = importlib.util.spec_from_file_location('historical_calibration_trainer', HERE / 'source/run_m0.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.HERE = output
    module.ROOT = inputs.parent
    module.BP = inputs
    module.main()


if __name__ == '__main__':
    main()
