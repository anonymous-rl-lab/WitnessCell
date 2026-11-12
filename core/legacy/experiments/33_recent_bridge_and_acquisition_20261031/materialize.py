"""Restore every original experiment file, validating its original SHA256.

python materialize.py --out /absolute/new_directory
Uses only Python's standard library. No network, float conversion or pickle loading.
"""
import argparse, hashlib, json, os, sys, tarfile, zlib
from pathlib import Path, PurePosixPath

R=Path(__file__).resolve().parent
def sha(b):return hashlib.sha256(b).hexdigest()
def restore(out):
    out=Path(out).resolve()
    assert not out.exists(), 'Use a new output directory; existing results are never overwritten'
    blobs={}
    with tarfile.open(R/'assets.tar.xz','r:xz') as t:
        index=json.loads(t.extractfile('index.json').read())
        assert index['format']=='witnesscell-lossless-v1'
        for m in t.getmembers():
            if m.name.startswith('blobs/'):
                assert m.isfile() and len(m.name)==70
                b=t.extractfile(m).read();h=m.name.split('/')[1]
                codec=index.get('blob_codecs',{}).get(h)
                if codec:
                    assert codec['codec']=='byte_shuffle'
                    off,size=codec['offset'],codec['itemsize'];body=b[off:];n=len(body)//size
                    original=bytearray(len(body))
                    for j in range(size):original[j::size]=body[j*n:(j+1)*n]
                    b=b[:off]+original
                assert sha(b)==h,'Corrupt blob '+h;blobs[h]=b
    files={}
    for h,recipe in index['object_recipes'].items():
        if recipe['type']=='raw':b=blobs[recipe['blob']]
        else:
            assert recipe['type']=='zip_segments';parts=[]
            for p in recipe['parts']:
                if 'raw' in p:piece=blobs[p['raw']]
                else:
                    co=zlib.compressobj(p['level'],zlib.DEFLATED,-15)
                    piece=co.compress(blobs[p['deflate']])+co.flush()
                    assert sha(piece)==p['sha256_compressed'],(
                        'Deflate byte mismatch. Use zlib '+index['zlib_runtime']+' for exact archival reconstruction.')
                parts.append(piece)
            b=b''.join(parts)
        assert len(b)==recipe['size'] and sha(b)==h,'File reconstruction mismatch '+h
        files[h]=b
    out.mkdir(parents=True)
    paths={}
    for f in index['members']:
        p=PurePosixPath(f['path']);assert not p.is_absolute() and '..' not in p.parts
        if str(p) in paths:assert paths[str(p)]==f['sha256'],'Conflicting archive members'
        else:
            dest=out/str(p);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(files[f['sha256']])
            if f['mode']:os.chmod(dest,f['mode'])
            paths[str(p)]=f['sha256']
    result={'status':'PASS','original_archives':len(index['sources']),
            'original_member_records':len(index['members']),'restored_paths':len(paths),
            'unique_file_contents':len(files),'every_stored_original_file_sha256_verified':True,
            'cache_files_to_regenerate':len(index.get('regenerated_files',[])),
            'regeneration_command':'Use reproduce.py --out <fresh_directory> to restore all three cache families and verify all original files',
            'numeric_precision_changed':False,'zlib_runtime':zlib.ZLIB_RUNTIME_VERSION}
    (out/'RESTORATION_RESULT.json').write_text(json.dumps(result,indent=2))
    (out/'ORIGINAL_SOURCES.json').write_text(json.dumps(index['sources'],indent=2))
    (out/'REGENERATED_CACHE_MANIFEST.json').write_text(json.dumps(index.get('regenerated_files',[]),indent=2))
    print(json.dumps(result,indent=2));return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True)
    restore(ap.parse_args().out)
