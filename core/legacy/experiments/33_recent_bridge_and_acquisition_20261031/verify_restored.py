#!/usr/bin/env python3
"""Verify every original archive member by exact size and SHA256, never tolerance."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile

ROOT = Path(__file__).resolve().parent


def verify(restored, payload=ROOT / 'assets.tar.xz'):
    restored = Path(restored).resolve()
    with tarfile.open(payload, 'r:xz') as archive:
        index = json.load(archive.extractfile('index.json'))
    if index['format'] != 'witnesscell-lossless-v1':
        raise RuntimeError('Unsupported payload format')
    records = index['members'] + index.get('regenerated_files', [])
    if len(records) != 6838:
        raise RuntimeError(f'Expected 6838 original member records, found {len(records)}')
    verified = {}; per_archive = Counter(); total_bytes = 0
    for record in records:
        relative = PurePosixPath(record['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise RuntimeError('Unsafe member path: ' + str(relative))
        path = restored / str(relative)
        expected = (record['size'], record['sha256'])
        if str(relative) in verified:
            if verified[str(relative)] != expected:
                raise RuntimeError('Conflicting duplicate member: ' + str(relative))
        else:
            if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(restored):
                raise RuntimeError('Missing or unsafe original file: ' + str(relative))
            size = path.stat().st_size
            with path.open('rb') as stream:
                digest = hashlib.file_digest(stream, 'sha256').hexdigest()
            if (size, digest) != expected:
                raise RuntimeError('Original size/SHA256 mismatch: ' + str(relative))
            verified[str(relative)] = expected
        per_archive[record['source']] += 1
        total_bytes += record['size']
    for source in index['sources']:
        if per_archive[source['name']] != source['files']:
            raise RuntimeError('Source member count mismatch: ' + source['name'])
    return {'status': 'PASS', 'original_archive_count': len(index['sources']),
            'original_member_records_verified': len(records),
            'unique_original_paths_verified': len(verified),
            'stored_member_records': len(index['members']),
            'regenerated_member_records': len(index.get('regenerated_files', [])),
            'bytes_across_original_member_records': total_bytes,
            'per_archive_member_records': dict(sorted(per_archive.items())),
            'every_original_size_and_sha256_exact': True,
            'numeric_tolerance_used': False,
            'scope': 'All original archive file members; generated external audit reports are not original members.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restored', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    report = args.report.resolve()
    if report.is_relative_to(args.restored.resolve()) or report.is_relative_to(ROOT):
        parser.error('Write verification reports outside both restored originals and the package')
    try:
        result = verify(args.restored)
    except Exception as error:
        result = {'status': 'FAIL', 'error': str(error)}
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(json.dumps(result, indent=2) + '\n')
        raise
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
