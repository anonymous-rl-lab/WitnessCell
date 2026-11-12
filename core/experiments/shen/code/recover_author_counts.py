#!/usr/bin/env python3
"""Map immutable author count files to the locked Shen interchange schema.

Only design annotations and raw count columns are used. No published fitness or
interaction scores enter this preparation. Source choices are made before fitting.
"""
from pathlib import Path
import csv, hashlib, json, re, datetime, zipfile
import xml.etree.ElementTree as ET
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'source'
COMMIT = '7c05a5a551dc5cd3a3c841d2f22215a450b353a3'
REPO = 'https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline'
PREFIX = 'notebooks/abundance_threshold_heuristic_development/mali_combined_counts/'
FILES = [
 '20160403_Hela_CV4_19mer_1mm_py_20160616222235_counts_combined.txt',
 '20160510_A549_19mer_1mm_py_20160616113153_counts_combined.txt',
 '20160627_HeLa_A549_CV4_19mer_1mm_py_20160628114018_counts_combined.txt',
 '20160803_293_CV4_19mer_1mm_py_20160805172522_counts_combined.txt',
]
SCHEMA = ['cell_line','replicate','day','construct_id','gene_a','guide_a','gene_b','guide_b','count']

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def alias(s): return s.replace('NonTargetingControlGuideForHuman','NonTargeting')
def write(p,x): p.write_text(json.dumps(x,indent=2)+'\n')

def design():
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(SRC/'Shen2017_TableS2.xlsx') as z:
        assert z.testzip() is None
        ss=[''.join(v.itertext()) for v in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        sheet=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
        ans={}; duplicate=0
        for row in sheet.findall('.//s:row',ns)[1:]:
            d={}
            for c in row.findall('s:c',ns):
                v=c.find('s:v',ns)
                if v is not None:
                    key=re.sub(r'\d','',c.attrib['r'])
                    d[key]=ss[int(v.text)] if c.attrib.get('t')=='s' else v.text
            if not d.get('K'): continue
            cid=alias(d['K']); probes=cid.split('__'); assert len(probes)==2
            genes=[alias(d[k]) for k in ('B','F')]
            genes=['NTC' if x.startswith('NonTargeting') else x for x in genes]
            anno=(genes[0],probes[0],genes[1],probes[1])
            if cid in ans:
                assert ans[cid]==anno, ('conflicting design',cid)
                duplicate+=1
            ans[cid]=anno
    return ans,duplicate

def main():
    out=ROOT/'measurement'
    if out.exists(): raise RuntimeError('Fresh measurement directory required')
    out.mkdir()
    annotations,duplicates=design()
    receipts={Path(x['file']).name:x for x in json.loads((SRC/'historical_counts_receipts.json').read_text())}
    sources=[];mapping=[];excluded=[];counts=Counter();genes=set();pairids=set();sampleids={};seen=set()
    output=out/'canonical_counts.csv'
    with output.open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(SCHEMA)
        for index,name in enumerate(FILES):
            path=SRC/'historical_counts'/name
            assert sha(path)==receipts[name]['sha256']
            sources.append(dict(file=str(path.relative_to(ROOT)),sha256=sha(path),
                git_blob=receipts[name]['git_blob'],commit=COMMIT,
                url=REPO+'/blob/'+COMMIT+'/'+PREFIX+name,
                download_url='https://raw.githubusercontent.com/ucsd-ccbb/mali-dual-crispr-pipeline/'+COMMIT+'/'+PREFIX+name))
            with path.open() as f:
                rd=csv.DictReader(f,delimiter='\t'); cols=[]
                for col in rd.fieldnames[1:]:
                    m=re.search(r'[_-]d(\d+)[_-](\d+)_',col); assert m,col
                    day,rawrep=map(int,m.groups())
                    cell='HeLa' if col.lower().startswith('hela') else ('A549' if col.startswith('A549') else '293T')
                    # The author scoring input identifies the later A549 d28 run.
                    # Keep its d3/d14 columns from the first run, d21/d28 from June.
                    if index==1 and day==28:
                        excluded.append({'file':name,'column':col,'reason':'Earlier d28 run is not the d28 column used in author A549 scoring input'})
                        continue
                    rep=({3:1,4:2}[rawrep] if cell=='293T' else rawrep)
                    assert rep in (1,2) and day in (3,14,21,28)
                    cols.append((col,cell,rep,day));mapping.append(dict(file=name,column=col,cell_line=cell,source_repeat=rawrep,replicate=rep,day=day))
                for r in rd:
                    cid=alias(r['construct_id']); assert cid in annotations,('unknown',cid)
                    ga,pa,gb,pb=annotations[cid]
                    genes.update(x for x in (ga,gb) if x!='NTC')
                    if ga!=gb and 'NTC' not in (ga,gb):pairids.add(tuple(sorted((ga,gb))))
                    for col,cell,rep,day in cols:
                        key=(cell,rep,day,cid);assert key not in seen;seen.add(key)
                        val=float(r[col]); assert val>=0 and val.is_integer()
                        writer.writerow([cell,rep,day,cid,ga,pa,gb,pb,int(val)])
                        counts[(cell,rep,day)]+=1;sampleids.setdefault((cell,rep,day),set()).add(cid)
    assert len(counts)==24 and len(genes)==73 and len(pairids)==2628
    assert set(counts.values())=={24908}
    assert len(set(frozenset(v) for v in sampleids.values()))==1
    # Independently anchor all eight selected A549 columns to the author scoring input.
    known=SRC/'recovered_counts/A549_CV4_counts_w_everything.txt'
    with known.open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    k={alias(r['construct_id']):r for r in rows};align=[]
    for x in mapping:
        if x['cell_line']!='A549':continue
        with (SRC/'historical_counts'/x['file']).open() as f:
            rr=list(csv.DictReader(f,delimiter='\t'))
        target=[c for c in rows[0] if re.search(r'[_-]d'+str(x['day'])+r'[_-]'+str(x['replicate'])+r'_',c)]
        assert len(target)==1
        matched=[r for r in rr if alias(r['construct_id']) in k]
        equal=sum(float(r[x['column']])==float(k[alias(r['construct_id'])][target[0]]) for r in matched)
        assert equal==len(matched)==24309
        align.append(dict(source_column=x['column'],author_column=target[0],matched_constructs=equal))
    ledger=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_commit=COMMIT,source_files=sources,source_sha256=hashlib.sha256(json.dumps(sources,sort_keys=True).encode()).hexdigest(),
        source_name='Shen2017 CV4 raw counts recovered from the official author repository history',
        source_url=REPO+'/tree/'+COMMIT,
        sheet_columns={'design':'S2A: B=gene A, F=gene B, K=SequenceID; raw count mapping below','count_columns':mapping},
        repeat_mapping={'HeLa':'source 1/2 retained','A549':'source 1/2 retained','293T':'source 293_CV4 replicate identifiers 3/4 relabelled 1/2; distinct originals retained'},
        control_identity_mapping='Explicit NonTargeting* gene annotations in S2A; author NonTargetingControlGuideForHuman alias maps to NonTargeting',
        units='control-relative log2 growth slope per day',declared_raw_counts=True,
        excluded_columns=excluded,design_duplicate_rows=duplicates,
        counts_sha256=sha(output),design_sha256=sha(SRC/'Shen2017_TableS2.xlsx'),
        raw_rows=sum(counts.values()),genes=len(genes),design_gene_pairs=len(pairids),
        samples=[dict(cell_line=c,replicate=r,day=d,constructs=n) for (c,r,d),n in sorted(counts.items())],
        A549_author_input_alignment=align,
        source_identity_notes=['HeLa CV4 is used; Trex2, CV4R4R, CRV4 and MV4 datasets are excluded by assay identity.',
          'Historical 293_CV4 is mapped to the 293T context reported for this paper; original filenames and sample labels are retained.',
          'S2A contains 600 dual-NTC design records; the count files contain 599 dual-NTC constructs.',
          'Public pooled S2B scores are not model inputs or query outcomes.'],
        real_model_fits_so_far=0,formal_query_scoring_so_far=0,
        protocol_amendment=False,preflight_thresholds_unchanged=True)
    write(out/'SOURCE_PROVENANCE.json',ledger)
    print(json.dumps({k:ledger[k] for k in ['raw_rows','genes','design_gene_pairs','counts_sha256','real_model_fits_so_far']},indent=2))
if __name__=='__main__':main()
