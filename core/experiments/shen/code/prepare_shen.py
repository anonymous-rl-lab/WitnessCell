#!/usr/bin/env python3
"""Prepare canonical count data WITHOUT training or scoring formal queries.

The actual publisher workbook mapping must be audited separately. This script
never infers which workbook columns are repeats, fitness, or GI scores.
"""
from pathlib import Path
import argparse,json,hashlib,datetime,itertools
import numpy as np
from shen_adapter import read_counts,construct_slopes,aggregate,split_identities,make_input,projection_info,NTC,CELLS

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,obj):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,indent=2)+'\n')
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts',type=Path,required=True,help='canonical source counts CSV; not the publisher XLSX')
    p.add_argument('--provenance',type=Path,required=True,help='audited field mapping and source hash JSON')
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    if args.out.exists():p.error('fresh output path required')
    provenance=json.loads(args.provenance.read_text())
    required={'source_sha256','source_name','source_url','sheet_columns','repeat_mapping','control_identity_mapping','units','declared_raw_counts'}
    if not required<=set(provenance) or provenance['declared_raw_counts'] is not True:
        p.error('Missing audited raw-count provenance; do not map published GI/whole-map singles')
    if not len(str(provenance['source_sha256']))==64:p.error('Source SHA256 missing')
    rows=read_counts(args.counts)
    cells=set(r['cell_line'] for r in rows)
    if cells!=set(CELLS):p.error('Keep all three source contexts; partial source must have a protocol amendment')
    identities={tuple(sorted((r['gene_a'],r['gene_b']))) for r in rows if r['gene_a']!=r['gene_b'] and NTC not in (r['gene_a'],r['gene_b'])}
    split=split_identities(identities) # source identities, BEFORE outcome-dependent coverage filters
    slopes,qc=construct_slopes(rows);tables,ex=aggregate(slopes)
    args.out.mkdir()
    write(args.out/'IDENTITY_SPLIT.json',split)
    write(args.out/'PREPROCESSING_QC.json',qc)
    write(args.out/'COVERAGE_EXCLUSIONS.json',ex)
    summary={}
    for cell in CELLS:
        if cell not in tables:p.error(f'No qualified two-repeat input for {cell}; no substitutions')
        table=tables[cell]
        inp=make_input(table,split['train'],split['query'])
        write(args.out/'training_inputs'/f'{cell}.json',inp)
        truth=np.array([table['pairs'][tuple(pair)] for pair in inp['queries']])
        dst=args.out/'sealed_evaluation'/f'{cell}.npz';dst.parent.mkdir(exist_ok=True)
        np.savez_compressed(dst,genes=np.array(inp['genes']),pairs=np.array(inp['queries']),y=truth)
        _,_,rank=projection_info(inp['genes'],inp['queries'])
        summary[cell]=dict(genes=len(inp['genes']),train_pairs=len(inp['training']),queries=len(inp['queries']),
                           partner_dimensions=len(inp['queries'])-rank,sealed_truth_sha256=sha(dst))
        if cell=='HeLa':
            sm=make_input(table,split['smoke_train'],split['smoke_validation'])
            write(args.out/'smoke_input/HeLa.json',sm)
            np.savez_compressed(args.out/'smoke_input/HeLa_truth.npz',genes=np.array(sm['genes']),pairs=np.array(sm['queries']),
                                y=np.array([table['pairs'][tuple(pair)] for pair in sm['queries']]))
    write(args.out/'PREPARATION_RECEIPT.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
          counts_sha256=sha(args.counts),provenance=provenance,summary=summary,
          source_values_read_for_measurement_preprocessing=True,
          formal_model_fitted=False,formal_predictions_scored=False,
          note='Sealing separates software access, not a cryptographic or institutional blind'))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
