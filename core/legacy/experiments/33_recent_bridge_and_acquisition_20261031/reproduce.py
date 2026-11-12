#!/usr/bin/env python3
"""Restore all nine rounds, rebuild exact caches, and replay historical results.

Usage: python -B reproduce.py --out /absolute/new_directory
The output must be new and outside this package. No original result is overwritten.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import zlib

ROOT = Path(__file__).resolve().parent
THREAD_KEYS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
               'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists() or out.is_relative_to(ROOT):
        parser.error('--out must be a nonexistent directory outside the package')
    out.mkdir(parents=True)
    logs = out / 'logs'; logs.mkdir()
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env['PYTHONOPTIMIZE'] = '0'
    env.update({key: '1' for key in THREAD_KEYS})
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    commands = []; started = time.monotonic()
    result = {'status': 'RUNNING', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Reconstruct and replay the nine original rounds; no new tuning or results.',
              'commands': commands, 'output': str(out)}

    def save():
        result['elapsed_seconds'] = time.monotonic() - started
        (out / 'COMPLETE_REPRODUCTION.json').write_text(json.dumps(result, indent=2) + '\n')

    def run(name, script, *arguments, cwd=ROOT):
        log = logs / (name + '.log')
        command = [sys.executable, '-B', str(script), *map(str, arguments)]
        stamp = time.monotonic()
        with log.open('w') as stream:
            process = subprocess.run(command, cwd=cwd, env=env, stdout=stream,
                                     stderr=subprocess.STDOUT)
        entry = {'stage': name, 'command': command, 'returncode': process.returncode,
                 'seconds': time.monotonic() - stamp, 'log': str(log.relative_to(out))}
        commands.append(entry)
        print(json.dumps({'stage': name, 'returncode': process.returncode}), flush=True)
        if process.returncode:
            raise RuntimeError(f'{name} failed (exit {process.returncode}); see {log}\n' +
                               log.read_text(errors='replace')[-4000:])

    def read(name):
        data = json.loads((out / name).read_text())
        if data.get('status') != 'PASS':
            raise RuntimeError('Stage report did not pass: ' + name)
        return data

    try:
        pins = json.loads((ROOT / 'ENVIRONMENT.json').read_text())
        actual = {'python': platform.python_version(),
                  'numpy': importlib.metadata.version('numpy'),
                  'scipy': importlib.metadata.version('scipy'),
                  'zlib_runtime': zlib.ZLIB_RUNTIME_VERSION, 'platform': platform.platform(),
                  'blas_threads': 1, 'bytecode_writes': False}
        expected = {'python': pins['python_required_for_verified_archival_bytes'],
                    'numpy': pins['numpy'], 'scipy': pins['scipy'],
                    'zlib_runtime': pins['zlib_runtime']}
        mismatches = {key: {'expected': value, 'actual': actual[key]}
                      for key, value in expected.items() if value != actual[key]}
        environment = {'actual': actual, 'required_exact_versions': expected,
                       'mismatches': mismatches}
        (out / 'ENVIRONMENT_CHECK.json').write_text(json.dumps(environment, indent=2) + '\n')
        result['environment_report'] = 'ENVIRONMENT_CHECK.json'
        if mismatches:
            raise RuntimeError('Exact reproduction environment differs from pins; see ENVIRONMENT_CHECK.json. '
                               'Install the stated Python, NumPy, SciPy and zlib versions; never edit frozen hashes.')
        restored = out / 'restored'
        run('01_materialize', ROOT / 'materialize.py', '--out', restored)
        bridge = restored / 'WitnessCell_NM_NCS_Bridge_20261031'
        jobs = [
            ('02_inner_cache', 'rebuild_inner_cache.py', ['--root', bridge, '--restore', '--report', out / 'INNER.json']),
            ('03_synthetic_cache', 'rebuild_synthetic_cache.py', ['--restored', restored, '--restore', '--report', out / 'SYNTHETIC.json']),
            ('04_jurkat_cache', 'rebuild_jurkat_cache.py', ['--bridge', bridge, '--restore', '--report', out / 'JURKAT.json'])]
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = [executor.submit(run, name, ROOT / script, *arguments)
                       for name, script, arguments in jobs]
            errors = []
            for future in futures:
                try:
                    future.result()
                except Exception as error:
                    errors.append(str(error))
            if errors:
                raise RuntimeError('\n'.join(errors))
        for report in ('INNER.json', 'SYNTHETIC.json', 'JURKAT.json'):
            read(report)
        run('05_exact_original_files', ROOT / 'verify_restored.py', '--restored', restored,
            '--report', out / 'EXACT_FILE_VERIFICATION.json')
        before = read('EXACT_FILE_VERIFICATION.json')
        run('06_replay_nine_rounds', ROOT / 'replay_suite.py', '--restored', restored,
            '--out', out / 'replay')
        suite = read('replay/REPLAY_SUITE_RESULT.json')
        replay_bridge = out / 'replay' / 'bridge_experiment_20261031'
        run('07_original_verify_release', replay_bridge / 'verify_release.py', cwd=replay_bridge)
        run('08_original_files_still_exact', ROOT / 'verify_restored.py', '--restored', restored,
            '--report', out / 'EXACT_FILE_VERIFICATION_AFTER_REPLAY.json')
        after = read('EXACT_FILE_VERIFICATION_AFTER_REPLAY.json')
        if before != after:
            raise RuntimeError('Original verification summary changed after replay')
        result.update(status='PASS', original_member_records_verified=before['original_member_records_verified'],
                      unique_original_paths_verified=before['unique_original_paths_verified'],
                      all_original_files_size_and_sha256_exact=True,
                      original_files_unchanged_after_replay=True,
                      replayed_stages=suite['requested_stages'],
                      original_release_verifier='PASS',
                      reports=['INNER.json', 'SYNTHETIC.json', 'JURKAT.json',
                               'EXACT_FILE_VERIFICATION.json', 'replay/REPLAY_SUITE_RESULT.json',
                               'EXACT_FILE_VERIFICATION_AFTER_REPLAY.json'],
                      archival_comparison='Exact original file size and SHA256; no numeric tolerance',
                      score_comparison='rtol=1e-12, atol=1e-12 where numerical score replay requires tolerance; see replay reports')
        save()
        print(json.dumps({'status': 'PASS', 'report': str(out / 'COMPLETE_REPRODUCTION.json'),
                          'original_member_records_verified': before['original_member_records_verified']}, indent=2))
    except Exception as error:
        result.update(status='FAIL', error=str(error))
        save()
        raise


if __name__ == '__main__':
    main()
