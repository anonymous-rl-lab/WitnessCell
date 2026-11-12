"""Lossless file/array deduplication. Reconstructs every original member byte-for-byte."""
from pathlib import Path, PurePosixPath
import argparse, ast, base64, hashlib, io, json, lzma, re, struct, tarfile, zipfile, zlib, time

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--sources',type=Path,required=True,help='Directory containing the nine original WitnessCell ZIP files')
ap.add_argument('--out',type=Path,required=True)
args=ap.parse_args();SRC=args.sources.resolve();DEST=args.out.resolve()
DEST.mkdir(parents=True,exist_ok=True)
blobs={}; recipes={}; members=[]; sources=[]; regenerated=[]; counters={'npz':0,'deflate_streams':0,'deflate_exact':0}
def sha(b):return hashlib.sha256(b).hexdigest()
def raw(b):
    h=sha(b);blobs.setdefault(h,b);return h
def encode(b,is_npz):
    h=sha(b)
    if h in recipes:return h
    if not is_npz:
        recipes[h]={'type':'raw','blob':raw(b),'size':len(b)};return h
    counters['npz']+=1
    try:
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            infos=sorted(z.infolist(),key=lambda x:x.header_offset)
            assert infos and all(not (i.flag_bits&1) for i in infos)
            parts=[];pos=0
            for i in infos:
                off=i.header_offset
                assert b[off:off+4]==b'PK\x03\x04'
                fn,extra=struct.unpack_from('<HH',b,off+26)
                begin=off+30+fn+extra;end=begin+i.compress_size
                assert pos<=begin<=end<=len(b)
                parts.append({'raw':raw(b[pos:begin])})
                compressed=b[begin:end]
                if i.compress_type==zipfile.ZIP_DEFLATED:
                    counters['deflate_streams']+=1
                    data=z.read(i);co=zlib.compressobj(6,zlib.DEFLATED,-15)
                    rebuilt=co.compress(data)+co.flush()
                    if rebuilt==compressed:
                        parts.append({'deflate':raw(data),'level':6,'sha256_compressed':sha(compressed)})
                        counters['deflate_exact']+=1
                    else:parts.append({'raw':raw(compressed)})
                else:parts.append({'raw':raw(compressed)})
                pos=end
            parts.append({'raw':raw(b[pos:])})
            recipes[h]={'type':'zip_segments','parts':parts,'size':len(b)}
    except (ValueError,AssertionError,zipfile.BadZipFile,NotImplementedError):
        recipes[h]={'type':'raw','blob':raw(b),'size':len(b)}
    return h

for p in sorted(SRC.glob('WitnessCell_*.zip')):
    source={'name':p.name,'sha256':sha(p.read_bytes()),'size_bytes':p.stat().st_size,'files':0,'expanded_bytes':0}
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        for i in z.infolist():
            if i.is_dir():continue
            q=PurePosixPath(i.filename)
            assert not q.is_absolute() and '..' not in q.parts
            b=z.read(i)
            training_cache=('/growth_run/training_v2/' in i.filename and (
                '/inner_models/' in i.filename or
                (i.filename.endswith('_inner_predictions.npz') and
                 not i.filename.endswith('/rank1_t0.3_n0.1_inner_predictions.npz')) or
                i.filename.endswith('/predictions.npz')))
            synthetic_cache=(i.filename=='WitnessCell_Formal_Hidden_Mechanism_Bridge_20261031/frozen/frozen_arrays.npz' or
                ('/causal_run/formal/' in i.filename and q.name in ('weights.npz','predictions.npz')))
            jurkat_cache=('/growth_assets/extracted/WitnessCell_Horlbeck_Jurkat_Confirmation_20261018/results/high/' in i.filename and
                (('/inner_models/' in i.filename and q.suffix=='.npz') or i.filename.endswith('_inner_predictions.npz')))
            if training_cache or synthetic_cache or jurkat_cache:
                regenerated.append({'source':p.name,'path':i.filename,'sha256':sha(b),
                                    'size':len(b),'rebuild':('rebuild_inner_cache.py' if training_cache else
                                        'rebuild_synthetic_cache.py' if synthetic_cache else 'rebuild_jurkat_cache.py')})
                source['files']+=1;source['expanded_bytes']+=len(b)
                continue
            h=encode(b,q.suffix=='.npz')
            members.append({'source':p.name,'path':i.filename,'sha256':h,'size':len(b),
                            'mode':(i.external_attr>>16)&0o777})
            source['files']+=1;source['expanded_bytes']+=len(b)
    sources.append(source)
    print(json.dumps({'source_encoded':p.name,'blob_count':len(blobs),'unique_bytes':sum(map(len,blobs.values()))}),flush=True)

codecs={}
def shuffled(h,b):
    if not b.startswith(b'\x93NUMPY'):return b
    try:
        version=b[6]
        start=10 if version==1 else 12
        n=struct.unpack_from('<H' if version==1 else '<I',b,8)[0]
        meta=ast.literal_eval(b[start:start+n].decode('latin1').strip())
        dtype=meta['descr'];m=re.fullmatch(r'[<>=|]?([fiuc])(\d+)',dtype)
        if not m:return b
        size=int(m.group(2));offset=start+n;body=b[offset:]
        if size<2 or size>16 or len(body)%size:return b
        codecs[h]={'codec':'byte_shuffle','offset':offset,'itemsize':size}
        return b[:offset]+b''.join(body[j::size] for j in range(size))
    except (ValueError,TypeError,SyntaxError,KeyError,struct.error):return b
packed_blobs={h:shuffled(h,b) for h,b in blobs.items()}
index={'format':'witnesscell-lossless-v1','compression':'xz solid + exact file/NPY deduplication + reversible byte shuffle',
       'numeric_precision_changed':False,'discarded_without_rebuild':False,'sources':sources,'members':members,
       'regenerated_files':regenerated,
       'object_recipes':recipes,'blob_codecs':codecs,'zlib_runtime':zlib.ZLIB_RUNTIME_VERSION}
(DEST/'CATALOG.json').write_text(json.dumps({'sources':sources,'total_original_files':len(members)+len(regenerated),
        'total_original_bytes':sum(x['expanded_bytes'] for x in sources),'unique_original_files':len(recipes),
        'regenerated_cache_files':len(regenerated),'regenerated_cache_bytes':sum(x['size'] for x in regenerated),
        'blob_count':len(blobs),'blob_uncompressed_bytes':sum(map(len,blobs.values())),**counters},indent=2))
(DEST/'REGENERATED_CACHE_MANIFEST.json').write_text(json.dumps(regenerated,indent=2))
payload=DEST/'assets.tar.xz.tmp'
with lzma.open(payload,'wb',preset=9) as stream:
    with tarfile.open(fileobj=stream,mode='w|',format=tarfile.PAX_FORMAT) as t:
        def add(name,b):
            inf=tarfile.TarInfo(name);inf.size=len(b);inf.mode=0o644;inf.mtime=0;t.addfile(inf,io.BytesIO(b))
        add('index.json',json.dumps(index,separators=(',',':'),sort_keys=True).encode())
        # Adjacent array-like blobs improve solid compression without altering bytes.
        for h,b in sorted(packed_blobs.items(),key=lambda x:(x[1][:16],len(x[1]),x[0])):add('blobs/'+h,b)
payload.replace(DEST/'assets.tar.xz');payload=DEST/'assets.tar.xz'
print(json.dumps({'payload':str(payload),'bytes':payload.stat().st_size,**counters}),flush=True)
