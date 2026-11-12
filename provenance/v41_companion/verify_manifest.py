"""Verify delivered supplement files before replay."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for row in manifest['files']:
    p=root/row['path']
    assert p.resolve().is_relative_to(root),row['path']
    assert p.is_file() and p.stat().st_size==row['bytes'],row['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
print('PASS',len(manifest['files']),'delivered files')
