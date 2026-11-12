"""Verify split-file checksums and atomically reconstruct the source H5AD.
Only checks transfer integrity; it does not authenticate the GEO original.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--parts-dir', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    a = p.parse_args()
    out = a.root.resolve(); (out/'audit').mkdir(parents=True, exist_ok=True)
    manifest_bytes = a.manifest.read_bytes(); m = json.loads(manifest_bytes)
    if Path(m['source_name']).name != m['source_name']:
        raise ValueError('Unsafe source filename')
    audit = {'source_manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(), 'parts': []}
    h, md = hashlib.sha256(), hashlib.md5()
    offset = 0; target = out / m['source_name']; tmp = target.with_suffix(target.suffix+'.assembling')
    try:
        with tmp.open('wb') as dest:
            for i, part in enumerate(sorted(m['parts'], key=lambda x: x['index']), 1):
                if part['index'] != i or part['offset'] != offset:
                    raise ValueError('Missing, duplicated or noncontiguous part')
                name = part['name']
                if Path(name).name != name:
                    raise ValueError('Unsafe part filename')
                src = a.parts_dir/name
                if not src.exists():
                    src = a.parts_dir/(name+'.bin')
                sh, mh, size = hashlib.sha256(), hashlib.md5(), 0
                with src.open('rb') as f:
                    for block in iter(lambda: f.read(8 << 20), b''):
                        sh.update(block);mh.update(block);h.update(block);md.update(block)
                        dest.write(block);size += len(block)
                if (size != part['size_bytes'] or sh.hexdigest() != part['sha256']
                        or mh.hexdigest() != part['md5']):
                    raise ValueError(f'Checksum/size mismatch: {name}')
                audit['parts'].append({'name':name,'size':size,'sha256':sh.hexdigest(),'md5':mh.hexdigest(),'pass':True})
                offset += size
        if offset != m['size_bytes'] or h.hexdigest() != m['sha256'] or md.hexdigest() != m['md5']:
            raise ValueError('Reconstructed source checksum/size mismatch')
        os.replace(tmp, target)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    audit.update(size_bytes=offset,sha256=h.hexdigest(),md5=md.hexdigest(),
                 status=f'PASS_ALL_{len(audit["parts"])}_PARTS_AND_REASSEMBLED_FILE',
                 reference='User transfer manifest and its original Drive metadata; not a GEO publisher checksum comparison')
    (out/'audit/TRANSFER_VERIFIED.json').write_text(json.dumps(audit,indent=2))
    print(json.dumps({k:v for k,v in audit.items() if k != 'parts'},indent=2))

if __name__ == '__main__':
    main()
