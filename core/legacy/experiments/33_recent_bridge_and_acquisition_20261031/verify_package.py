"""Verify the compact package itself before restoring historical source files."""
from pathlib import Path
import hashlib,json

root=Path(__file__).resolve().parent
manifest=json.loads((root/'PACKAGE_MANIFEST.json').read_text())
expected={r['path']:r for r in manifest['files']}
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()
        and '__pycache__' not in p.parts and p.name!='PACKAGE_MANIFEST.json'}
assert actual==set(expected),('File set mismatch', sorted(actual-set(expected)),sorted(set(expected)-actual))
for name,row in expected.items():
    p=root/name
    assert p.stat().st_size==row['bytes'],('Size mismatch',name)
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],('Hash mismatch',name)
print(json.dumps({'status':'PASS','compact_package_files':len(expected),'original_sources':9}))
