"""Archive each frozen job once; replay, checksums and durable receipts precede success.

No training, budget, split or model-selection changes. Never blindly retry an
upload with an unknown outcome. The index is updated only after confirmed writes.
"""
from pathlib import Path
import argparse, fcntl, hashlib, json, os, subprocess, sys, time, zipfile

R = Path(__file__).resolve().parent
D = R / 'completed_archives_20261106'
INDEX = D / 'WitnessCell_Completed_Batches_20261106.json'
IDENTITY = D / 'INDEX_ID.json'
CORE = ('train_cell_models.py', 'train_cell_models_runtime.py', 'lazy_cell_data.py', 'fast_gears_ops.py')
BASE = 'libfile_88e02f619d7c8191af9724c6000ff85a'
CELLS = 'libfile_1d528a9a0cb88191ac870f34a63e56a6'

def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + '.tmp')
    with tmp.open('w') as f:
        json.dump(value, f, indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
    tmp.replace(p)

def read(p, default=None):
    return json.loads(p.read_text()) if p.exists() else default

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        while b := f.read(8 * 1024 * 1024): h.update(b)
    return h.hexdigest()

def frozen_jobs(root=R):
    for fold in range(5):
        for kind in ('gears', 'scouter'):
            for seed in range(3):
                p = root / 'results/detected49' / f'{kind}_seed{seed}' / f'fold_{fold}'
                if (p / 'FREEZE.json').exists():
                    yield dict(name=f'{kind}_seed{seed}_fold{fold}', kind=kind, seed=seed, fold=fold), p

def validate_frozen(p, root=R):
    f = read(p / 'FREEZE.json')
    assert f['status'] == 'FROZEN_UNSCORED', p
    assert f['outer_scores_read'] is False, p
    assert not set(f['train']).intersection(f['query']), p
    assert sha(p / 'predictions.npz') == f['predictions_sha256'], p
    assert read(p / 'refit/RUN.json')['status'] == 'COMPLETE', p
    assert (p / 'refit/last.pt').stat().st_size > 0, p
    for name, digest in f.get('runtime_fingerprint', {}).items():
        assert sha(root / name) == digest, f'Frozen source changed: {name}'
    return f

def verify_one(p):
    # Executed in a short-lived subprocess so native-model memory is released.
    from train_cell_models import torch, GM, SM, predict, load, np
    from lazy_cell_data import LazyData
    torch.set_num_threads(1)
    kind = p.parent.name.split('_seed')[0]
    f = validate_frozen(p)
    d = LazyData(f['fold'], kind); model = d.newmodel()
    assert type(model) is (GM.GEARS_Model if kind == 'gears' else SM.ScouterModel)
    model.load_state_dict(torch.load(p / 'refit/last.pt', weights_only=True)); model.eval()
    replay = predict(model, d, d.q); stored = load(p / 'predictions.npz')
    assert np.array_equal(stored['query'], d.q)
    assert np.array_equal(stored['genes'], d.m['genes'])
    np.testing.assert_allclose(replay, stored['pred_effect'], rtol=1e-4, atol=5e-5)
    write(p / 'BATCH_NATIVE_REPLAY.json', dict(status='PASS', original_upstream_forward=True,
        maximum_absolute_difference=float(np.max(abs(replay - stored['pred_effect']))),
        weight_sha256=sha(p / 'refit/last.pt'), prediction_sha256=f['predictions_sha256'],
        outer_outcomes_scored=False, checked_unix=time.time()))

def build(job, p):
    f = validate_frozen(p); weight = sha(p / 'refit/last.pt')
    proof = None
    for q in sorted(p.glob('*REPLAY*.json')):
        a = read(q)
        if a.get('status') == 'PASS' and a.get('weight_sha256') == weight and (a.get('original_upstream_forward') or a.get('native_source')):
            if a.get('prediction_sha256', f['predictions_sha256']) == f['predictions_sha256']:
                proof = a; break
    if proof is None:
        subprocess.run([sys.executable, str(Path(__file__)), '--verify-one', str(p)], cwd=R, check=True)
        proof = read(p / 'BATCH_NATIVE_REPLAY.json')
    files = sorted({q for q in p.rglob('*') if q.is_file() and (q.suffix == '.json' or q.name == 'predictions.npz' or q == p / 'refit/last.pt')})
    files += [R / name for name in CORE] + [R / 'PROTOCOL.md']
    dest = D / f"WitnessCell_Completed_{job['name']}_20261106.zip"
    manifest = {}
    with zipfile.ZipFile(dest.with_suffix('.tmp'), 'w', zipfile.ZIP_DEFLATED, compresslevel=1) as z:
        for q in files:
            name = str(q.relative_to(R)); z.write(q, name)
            manifest[name] = dict(bytes=q.stat().st_size, sha256=sha(q))
        z.writestr('MANIFEST.json', json.dumps(manifest, indent=2))
        z.writestr('RECOVERY.json', json.dumps(dict(job=job, base_code_library_file_id=BASE,
            cell_data_library_file_id=CELLS, instructions='Overlay on the base run directory. Apply completed job archives AFTER the latest partial checkpoint. Validate FREEZE, prediction hash and weight hash before scheduling. Completed jobs must be skipped.',
            outer_outcomes_scored=False), indent=2))
    dest.with_suffix('.tmp').replace(dest)
    # Independent archive read catches truncation and validates exact restored bytes.
    with zipfile.ZipFile(dest) as z:
        for name, info in manifest.items():
            h = hashlib.sha256(); size = 0
            with z.open(name) as stream:
                while b := stream.read(8 * 1024 * 1024): h.update(b); size += len(b)
            assert size == info['bytes'] and h.hexdigest() == info['sha256'], name
    print(json.dumps(dict(event='JOB_VERIFIED_AND_PACKED', job=job['name'], bytes=dest.stat().st_size)), flush=True)
    return dict(**job, relative_folder=str(p.relative_to(R)), freeze_sha256=sha(p / 'FREEZE.json'),
        prediction_sha256=f['predictions_sha256'], weight_sha256=weight, native_replay=proof,
        archive_sha256=sha(dest), archive_bytes=dest.stat().st_size, local_path=str(dest), verified_unix=time.time())

def upload(key, entries, metadata):
    """Durably journal one ordered mutation; confirmed receipts are restart-safe."""
    pending = D / f'{key}_PENDING.json'
    stdout = D / f'{key}_UPLOAD.stdout'
    stderr = D / f'{key}_UPLOAD.stderr'
    if pending.exists():
        saved = read(pending)
        receipt = parse_receipt(stdout)
        if receipt is None:
            raise RuntimeError(f'Upload outcome unknown; inspect {pending.name}; no automatic retry')
    else:
        saved = dict(uploads=entries, metadata=metadata, started_unix=time.time())
        write(pending, saved)
        helper = os.environ.get('WITNESSCELL_LIBRARY_UPLOAD_HELPER')
        if not helper or not Path(helper).is_file():
            pending.unlink(); raise RuntimeError('Current Library upload helper is required')
        with stdout.open('w') as out, stderr.open('w') as err:
            subprocess.run([sys.executable, helper], input=json.dumps(dict(uploads=entries)),
                text=True, stdout=out, stderr=err)
        receipt = parse_receipt(stdout)
        if receipt is None:
            raise RuntimeError(f'Upload outcome unknown; inspect {pending.name}; no automatic retry')
    rows = receipt['results']
    assert len(rows) == len(saved['uploads']), 'Receipt count mismatch'
    for expected, row in zip(saved['uploads'], rows):
        assert row.get('status') == 'succeeded', f'Upload not confirmed: {row}'
        assert row['local_path'] == expected['local_path'], 'Receipt path mismatch'
        assert row.get('local_metadata_applied') is True, 'Local identity not applied'
    return saved, rows

def parse_receipt(path):
    if not path.exists(): return None
    for line in reversed(path.read_text().splitlines()):
        try:
            obj = json.loads(line)
            if isinstance(obj, dict) and 'results' in obj: return obj
        except ValueError: pass
    return None

def empty_index():
    return dict(schema=1, total_jobs=30, jobs={}, base_code_library_file_id=BASE,
        cell_data_library_file_id=CELLS, instructions='Restore the base and latest in-progress snapshot first, then overlay these immutable completed-job archives. Check archive_sha256 and each MANIFEST.json. Valid completed jobs are excluded from training. Scores require all 30 jobs.')

def consume_jobs(index, saved, rows):
    for meta, row in zip(saved['metadata'], rows):
        index['jobs'][meta['name']] = dict(**meta, library_file_id=row['library_file_id'],
            file_id=row['file_id'], version=row['current_version_number'], archived_unix=time.time())
    index['updated_unix'] = time.time(); index['completed_archived'] = len(index['jobs'])
    write(INDEX, index)
    write(D / 'LAST_JOB_RECEIPT.json', dict(results=rows))
    (D / 'JOBS_PENDING.json').unlink()

def archive():
    D.mkdir(exist_ok=True)
    with (D / 'ARCHIVE.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('Archive worker already running.', flush=True); return
        index = read(INDEX, empty_index())
        if (D / 'JOBS_PENDING.json').exists():
            saved, rows = upload('JOBS', [], [])
            consume_jobs(index, saved, rows)
        if (D / 'INDEX_PENDING.json').exists():
            saved, rows = upload('INDEX', [], [])
            write(IDENTITY, dict(**rows[0], content_sha256=saved['metadata']['content_sha256']))
            (D / 'INDEX_PENDING.json').unlink()
        pending = []
        for job, p in frozen_jobs():
            prior = index['jobs'].get(job['name'])
            if prior:
                assert prior['freeze_sha256'] == sha(p / 'FREEZE.json'), 'Archived job freeze changed'
                continue
            pending.append((job, p))
        metadata = [build(job, p) for job, p in pending]
        if metadata:
            entries = [dict(local_path=m['local_path'], purpose='create_library_file', library_artifact_type='other') for m in metadata]
            saved, rows = upload('JOBS', entries, metadata)
            consume_jobs(index, saved, rows)
            print(json.dumps(dict(event='JOBS_ARCHIVED', jobs=[m['name'] for m in metadata], count=len(index['jobs']))), flush=True)
        if INDEX.exists():
            identity = read(IDENTITY, {})
            digest = sha(INDEX)
            if identity.get('content_sha256') != digest:
                entry = dict(local_path=str(INDEX), purpose='create_library_file', library_artifact_type='other')
                if identity:
                    entry = dict(local_path=str(INDEX), purpose='replace_library_file', library_file_id=identity['library_file_id'],
                        expected_current_version=identity['current_version_number'], version_reason='New verified completed experiment batch')
                saved, rows = upload('INDEX', [entry], dict(content_sha256=digest))
                write(IDENTITY, dict(**rows[0], content_sha256=digest))
                (D / 'INDEX_PENDING.json').unlink()
        write(D / 'STATUS.json', dict(status='ARCHIVED', completed_archived=len(index['jobs']),
            newly_archived=len(metadata), checked_unix=time.time(), index_identity=read(IDENTITY)))
        print(json.dumps(dict(status='ARCHIVED', completed_archived=len(index['jobs']), newly_archived=len(metadata))), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify-one', type=Path)
    args = parser.parse_args()
    if args.verify_one: verify_one(args.verify_one)
    else:
        try: archive()
        except Exception as exc:
            write(D / 'STATUS.json', dict(status='ARCHIVE_BLOCKED', error=str(exc), checked_unix=time.time()))
            raise
