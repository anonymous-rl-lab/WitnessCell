"""Verify and overlay downloaded completed-job archives, then list jobs to resume."""
from pathlib import Path, PurePosixPath
import argparse, hashlib, json, os, shutil, tempfile, zipfile
from archive_completed_jobs import sha, validate_frozen

def restore(index_path, archive_dir, destination, verify_only=False):
    index = json.loads(index_path.read_text()); destination = destination.resolve()
    restored = []
    for name, row in index['jobs'].items():
        archive = archive_dir / Path(row['local_path']).name
        assert sha(archive) == row['archive_sha256'], f'Archive checksum failed: {name}'
        with zipfile.ZipFile(archive) as z:
            manifest = json.loads(z.read('MANIFEST.json'))
            assert len(z.namelist()) == len(set(z.namelist())), 'Duplicate archive paths'
            assert set(z.namelist()) == set(manifest) | {'MANIFEST.json', 'RECOVERY.json'}
            for relative, info in manifest.items():
                path = PurePosixPath(relative)
                assert not path.is_absolute() and '..' not in path.parts
                target = destination / relative
                assert target.resolve().is_relative_to(destination), 'Unsafe archive path'
                h = hashlib.sha256(); size = 0
                with z.open(relative) as stream:
                    while block := stream.read(8 * 1024 * 1024): h.update(block); size += len(block)
                assert size == info['bytes'] and h.hexdigest() == info['sha256'], relative
            folder = destination / row['relative_folder']
            freeze = folder / 'FREEZE.json'
            was_frozen = freeze.exists()
            if was_frozen:
                assert sha(freeze) == row['freeze_sha256'], 'Existing frozen job differs; refusing overwrite'
            if not verify_only:
                for relative, info in manifest.items():
                    target = destination / relative
                    if target.exists() and sha(target) == info['sha256']: continue
                    if target.exists() and (relative.endswith('.py') or was_frozen):
                        raise RuntimeError(f'Conflicting existing file: {relative}')
                    target.parent.mkdir(parents=True, exist_ok=True)
                    # Each complete file replaces a partial file atomically.
                    with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as out:
                        tmp = Path(out.name)
                        with z.open(relative) as src: shutil.copyfileobj(src, out, 8 * 1024 * 1024)
                        out.flush(); os.fsync(out.fileno())
                    tmp.replace(target)
                validate_frozen(folder, root=destination)
                assert sha(folder / 'refit/last.pt') == row['weight_sha256']
            restored.append(name)
    remaining = [f'{k}_seed{s}_fold{f}' for f in range(5) for k in ('gears', 'scouter') for s in range(3)
        if f'{k}_seed{s}_fold{f}' not in restored]
    return dict(status='PASS', completed_skipped=restored, completed_count=len(restored),
        remaining_count=len(remaining), remaining=remaining, verification_only=verify_only)

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=Path,required=True)
    p.add_argument('--archive-dir',type=Path,required=True);p.add_argument('--destination',type=Path,required=True)
    p.add_argument('--verify-only',action='store_true');a=p.parse_args()
    print(json.dumps(restore(a.index,a.archive_dir,a.destination,a.verify_only),indent=2))
