#!/usr/bin/env python3
"""Assemble and check v24 SI using immutable copied source tables; no training."""
from pathlib import Path
import argparse, csv, hashlib, json, math, re, shutil, statistics

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'source_data'
P = SRC / 'publication'
F = SRC / 'formal_bridge'
G = SRC / 'gse220974'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def rows(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as h:
        return list(csv.DictReader(h))

def js(path):
    return json.loads(Path(path).read_text())

def table(headers, data):
    def cell(x):
        return str(x).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(map(cell, headers)) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(cell, r)) + ' |' for r in data])

def fmt(x, n=6):
    if x is None or str(x).strip() in ['', 'nan', 'NaN']:
        return '—'
    return f'{float(x):.{n}f}'

def snapshot(repo, bridge, main):
    records=[]
    def cp(src, dest, origin):
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
        records.append({'path':str(dest.relative_to(ROOT)), 'origin':origin,
                        'bytes':dest.stat().st_size, 'sha256':sha(dest)})
    for src in sorted((repo/'publication/source_data').rglob('*')):
        if src.is_file():
            cp(src, P/src.relative_to(repo/'publication/source_data'), str(src.relative_to(repo)))
    cp(repo/'publication/PUBLISHED_TABLES_v22.md', SRC/'legacy_tables_v22.md', 'publication/PUBLISHED_TABLES_v22.md')
    cp(repo/'REPRODUCTION_INDEX.json', SRC/'REPRODUCTION_INDEX.json', 'REPRODUCTION_INDEX.json')
    formal=bridge/'WitnessCell_Formal_Hidden_Mechanism_Bridge_20261031'
    names=['CONFIG.json','PROTOCOL.txt','DERIVED_REPORT_STATISTICS.json','BUDGET_AND_INTERVAL_REVIEW.json',
           'theory/THREE_GATE_DESIGN_ENUMERATION.json','theory/formal_derivation.txt',
           'theory/TARGET_RISK_LOWER_BOUND.txt','theory/MECHANISM_CLASS_BOUNDARY.txt']
    names += [str(x.relative_to(formal)) for x in (formal/'scores').glob('*.csv')]
    for name in names:
        cp(formal/name,F/name,'experiments/33_recent_bridge_and_acquisition_20261031/restored/'+formal.name+'/'+name)
    ex=repo/'experiments/34_gse220974_five_model_final'
    for name in ['SUMMARY.csv','MAIN_FIVE_MODELS.csv','CHOICES.csv','PAIRED_COMPARISONS.csv','GROUPS.json',
                 'FACTS.json','COMPLETION.json','EVALUATION.json','INDEPENDENT_SCORE_CHECK.json']:
        cp(ex/'reference'/name,G/name,str((ex/'reference'/name).relative_to(repo)))
    for name in ['PROTOCOL.md','PREPROCESS.json','PROTOCOL_FREEZE.json']:
        cp(ex/'protocol'/name,G/'protocol'/name,str((ex/'protocol'/name).relative_to(repo)))
    for name in ['COMPACTION_CONTRACT.json','ENVIRONMENT.json','EXTERNAL_ASSETS.json']:
        cp(ex/name,G/name,str((ex/name).relative_to(repo)))
    for kind in ['gears','scouter']:
        for seed in range(3):
            for fold in range(5):
                stem=Path(f'runs/{kind}_seed{seed}/fold_{fold}')
                for tail in ['SELECTION.json','refit/RUN.json']:
                    cp(ex/stem/tail,G/stem/tail,str((ex/stem/tail).relative_to(repo)))
    for stage, dirname, sub in [('development','target_aligned_run_20261031','calibration'),
                                 ('outer','target_aligned_outer_20261031','outer')]:
        for name in ['summary.json','rows.json','cases.json']:
            src=bridge/dirname/'scores'/sub/name
            cp(src,SRC/'acquisition'/stage/name,'experiments/33_recent_bridge_and_acquisition_20261031/restored/'+str(src.relative_to(bridge)))
    cp(main,ROOT/'WitnessCell_v24_Main_Methods.md','WitnessCell_v24_Main_Methods.md (unchanged companion)')
    (ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':records},indent=2)+'\n')

def legacy_blocks():
    text=(SRC/'legacy_tables_v22.md').read_text()
    parts=re.split(r'^### Supplementary Table S(\d+) \| ',text,flags=re.M)
    return {int(parts[i]):parts[i+1].strip() for i in range(1,len(parts),2)}

def clean_old(block):
    block=block.split('\n',1)[1].strip().replace('\\clearpage','')
    replacements={
        'Final v14':'Released core', 'final v14':'released core',
        'the separate matched RNA comparison is Table S2':'the separate matched RNA comparison is Table S10',
        'Note 2.3 explains the neural complete-response errors; Note 2.4 defines the distinct rectangle and projected outcomes.':'Notes 4 and 5 distinguish the response, rectangle and projected outcomes and describe the neural selection scope.',
        'see Note 2.3 for selection scope.':'see Note 5.2 for selection scope.',
        'Note 3.3':'Note 6.3',
        'Note 5.3':'Note 8.3',
        'Conditional sensitivity and adverse endpoint deletion are in Note 9.4.':'Conditional sensitivity and adverse endpoint deletion remain in the corresponding fixed-reference source records.',
    }
    for a,b in replacements.items():block=block.replace(a,b)
    return block.strip()

def build():
    old=legacy_blocks(); t={}
    rr=rows(F/'scores/summary.csv')
    labels={'OLS_low':'Candidate least squares, comparison panel','OLS_high':'Candidate least squares, informative panel',
            'contrast_GLS_high':'Contrast GLS, informative panel','zero_target':'Always predict zero'}
    t['FORMAL_RISK']=table(['Estimator','Mean target MSE','Monte Carlo SE','95% Monte Carlo interval'],
        [(labels[r['method']],fmt(r['four_world_rectangle_MSE']),fmt(r['MCSE']),f"[{fmt(r['MC95_low'])}, {fmt(r['MC95_high'])}]") for r in rr])
    t['FORMAL_RISK']+='\n\n'+table(['Paired comparison','Risk difference','95% Monte Carlo interval'],
        [(r['comparison'],fmt(r['paired_risk_difference']),f"[{fmt(r['MC95_low'])}, {fmt(r['MC95_high'])}]") for r in rows(F/'scores/paired_comparisons.csv')])
    wr=[r for r in rows(F/'scores/world_summary.csv') if r['method']=='OLS_high']
    t['FORMAL_WORLDS']=table(['World','Target MSE','Mean query MSE','Target coverage %','Empty set %','Nonzero declaration %'],
        [(r['world'],fmt(r['rectangle_MSE']),fmt(r['query_mean_MSE']),fmt(100*float(r['target_coverage']),1),
          fmt(100*float(r['empty_rate']),1),fmt(100*float(r['nonzero_certificate_rate']),1)) for r in wr])
    radius=math.sqrt(-2*math.log(.05))
    budget=[]
    for n in [64,128,256]:
        s=min(radius,max(math.sqrt(n/8)-radius,0))
        prob=1-math.exp(-s*s/2)
        budget.append((n,12*n,fmt(n/8,1),fmt(prob,5),'Yes' if prob>=.8 else 'No'))
    t['BUDGET']=table(['Observations per condition','Panel observations','Minimum squared line separation','Conservative probability','Meets 0.8 target'],budget)
    lb=rows(P/'released_benchmark/formal_top100_leaderboard.csv')
    t['LEADERBOARD']=table(['Position','Method','Mean operation rank','Mean metric rank','Records'],
        [(int(float(r['official_position'])),r['method'],fmt(r['mean_operation_rank']),fmt(r['mean_metric_rank']),r['operation_seed_units']) for r in lb])
    mm=rows(P/'released_benchmark/official_top100_per_metric_leaderboard.csv'); metrics=list(dict.fromkeys(r['metric'] for r in mm))
    lookup={(r['method'],r['metric']):r['position'] for r in mm}
    t['METRIC_POSITIONS']=table(['Method']+metrics,[(r['method'],*[lookup[r['method'],m] for m in metrics]) for r in lb])
    t['BENCHMARK_PAIRED']=table(['Comparator','WitnessCell minus comparator','95% cluster interval'],
        [(r['competitor'],fmt(r['rank_difference_witness_minus_competitor']),
          '['+', '.join(fmt(x) for x in json.loads(r['cluster_bootstrap_95_ci']))+']') for r in rows(P/'released_benchmark/paired_cluster_bootstrap.csv')])
    s7=clean_old(old[7]); t['BREADTH']=s7.split('**Separate analytic-risk')[0].strip()
    t['RISK_TABLES']='**Separate analytic-risk'+s7.split('**Separate analytic-risk',1)[1]
    t['RISK_TABLES']=t['RISK_TABLES'].split('### Supplementary')[0].strip()
    t['STRUCTURAL']=clean_old(old[1]);t['RNA_MATCHED']=clean_old(old[2]);t['RNA_CONTROLS']=clean_old(old[3])
    fold=rows(P/'Step2_ED5_Folds.csv'); fr=[]
    for f in range(5):
        ss={r['model']:r for r in fold if r['seed']=='0' and r['fold']==str(f)}
        b=ss['M0']; gain=lambda m,k:100*(1-float(ss[m][k])/float(b[k]))
        fr.append((f,b['conditions'],b['rectangles'],b['contrast_dimension'],fmt(gain('M1','C_SSE'),4),
                   fmt(gain('R','C_SSE'),4),fmt(gain('R','response_SSE'),4)))
    t['RNA_FOLDS']=table(['Fold','Pairs','Rectangles','Projected dimensions','M1 rectangle gain %','R rectangle gain %','R full gain %'],fr)
    t['RNA_RECTANGLES']=table(['Rectangle','Fold (zero-based)','Four pairs in signed order'],
       [(r['rectangle'],int(r['fold'])-1,'; '.join(r['ordered_pairs'].split('|'))) for r in rows(P/'RNA_rectangle_definition.csv')])
    t['GROWTH_K562']=clean_old(old[4]); t['GROWTH_JURKAT']=clean_old(old[5]);t['CHRONOLOGY']=clean_old(old[6])
    allg=rows(G/'SUMMARY.csv'); primary={'R_G','StandardRidge_RNA','SingleB_additive','GEARS_Ensemble','Scouter_Ensemble'}
    t['GSE_ALL']=table(['Archived model','Primary set','RNA MSE','Partner MSE','Mean regret','Hits / 28','Larger-pool regret'],
        [('PairPick (R-G) / `R_G`' if r['model']=='R_G' else '`'+r['model']+'`','Yes' if r['model'] in primary else '—',
          fmt(r['RNA_MSE']),fmt(r['partner_MSE']),fmt(r['mean_regret']),r['hits'],fmt(r['large_regret'])) for r in allg])
    t['GSE_PAIRED']=table(['Reference','R-G minus reference utility','95% anchor-bootstrap interval','Anchor clusters'],
        [(r['reference'],fmt(r['mean_paired_U_difference'],9),f"[{fmt(r['cluster_bootstrap_low'],9)}, {fmt(r['cluster_bootstrap_high'],9)}]",r['physical_anchor_clusters']) for r in rows(G/'PAIRED_COMPARISONS.csv')])
    gg=js(G/'GROUPS.json')
    t['GSE_GROUPS']=table(['Task','Fold','Mode-specific anchor','Candidates','Partners'],
        [(i+1,r['fold'],r['anchor'],r['n'],', '.join(r['partners'])) for i,r in enumerate(gg)])
    ep=[]
    for kind in ['gears','scouter']:
        for seed in range(3):
            epochs=[]; schedules=[]
            for f in range(5):
                d=js(G/f'runs/{kind}_seed{seed}/fold_{f}/SELECTION.json');s=d.get('selected',d)
                epochs.append(str(s['selected_epoch']));schedules.append(str(s.get('schedule','—')))
            ep.append((kind,seed,', '.join(epochs),', '.join(schedules) if kind=='gears' else 'Exponential factor 0.9'))
    t['GSE_EPOCHS']=table(['Model','Seed','Final-fit epochs, folds 0–4','LR halving interval, folds 0–4'],ep)
    ext=[]
    for label,ns in [('A. Frozen Dixit RNA and fold sensitivity',[13,14]),('B. Later Dixit program development',[15]),('C. Frozen Combi-seq comparison',[16]),('D. Combi-seq folds and experimental-run targets',[17])]:
        ext.append('**'+label+'**\n\n'+'\n\n'.join(clean_old(old[n]) for n in ns))
    t['EXTERNAL_TABLES']='\n\n'.join(ext)
    acq=[]
    for stage in ['development','outer']:
        d=js(SRC/'acquisition'/stage/'summary.json'); gains=d['gains']['calibrated']['uniform_expected']; each=d['gains_each']['calibrated']['uniform_expected']
        acq.append((stage,len(d['tasks']) if isinstance(d['tasks'],list) else d['tasks'],
                    fmt(100*gains['partner'],5),fmt(100*gains['AB'],5),fmt(100*gains['full'],5),
                    ', '.join(fmt(100*v,5) for v in each['partner'])))
    t['ACQUISITION']=table(['Stage','Targets','Partner gain %','Target-AB gain %','Complete gain %','Partner gains by repeat %'],acq)
    t['ACQUISITION']+='\n\nAll gains use the expected uniform acquisition reference and retain both measured repeats.'
    tasks=js(SRC/'REPRODUCTION_INDEX.json')['tasks']
    t['REPRODUCTION']=table(['Task name','Entry script','Recorded scope'],[(k,'`'+v['script']+'`',v['scope']) for k,v in tasks.items()])
    main=(ROOT/'WitnessCell_v24_Main_Methods.md').read_text()
    ref=main.split('## References',1)[1]
    wanted={3,4,6,7,10,12,13,14,15,16,17,18,19,20}
    t['REFERENCES']='\n\n'.join(re.sub(r'^(\d+)\.\s*',r'**[\1]** ',line) for line in ref.strip().splitlines() if (m:=re.match(r'^(\d+)\.',line)) and int(m.group(1)) in wanted)
    template=(ROOT/'authoring/SI_template.md').read_text()
    needed=set(re.findall(r'\{\{([A-Z_0-9]+)\}\}',template));assert needed==set(t),(needed-set(t),set(t)-needed)
    for k,v in t.items():template=template.replace('{{'+k+'}}',v)
    (ROOT/'WitnessCell_v24_Supplementary_Information.md').write_text(template)
    extra='# Additional archived development tables\n\nThese tables preserve v22 numbering and wording as source documentation, not v24 SI numbering. They are not additional independent confirmations. Original note references refer to the v22 supplement in the code archive.\n\n'
    extra+='\n\n'.join('## Archived Table S'+str(n)+' | '+old[n].replace('\\clearpage','') for n in [11,12,18,19])
    (ROOT/'additional_development_tables.md').write_text(extra)
    return template

def verify(text):
    import numpy as np
    from scipy.stats import t as student_t
    report={}
    manifest=js(ROOT/'SOURCE_MANIFEST.json')['files']
    bad=[r['path'] for r in manifest if sha(ROOT/r['path'])!=r['sha256']]
    assert not bad,bad
    report['source_files_checked']=len(manifest);report['source_hashes']='PASS'
    if (ROOT/'CLAIM_SOURCE_MAP.csv').exists():
        claims=rows(ROOT/'CLAIM_SOURCE_MAP.csv')
        assert len(claims)==35 and all((ROOT/r['source_path']).is_file() for r in claims)
        report['claim_source_rows_checked']=len(claims)
    ns=[int(x) for x in re.findall(r'^### Supplementary Table S(\d+) \|',text,re.M)]
    assert ns==list(range(1,25)),ns
    assert len(re.findall(r'^## Supplementary Note \d+ \|',text,re.M))==12
    assert '{{' not in text and '\\clearpage' not in text
    assert all(int(x)<=12 for x in re.findall(r'Note (\d+)\.',text))
    # Markdown table syntax: escaped pipes must not become extra columns.
    for block in re.findall(r'(?m)(?:^\|.*\n?)+',text):
        widths=[len(re.split(r'(?<!\\)\|',line))-2 for line in block.splitlines()]
        assert len(set(widths))==1,widths
    report['tables']=24;report['notes']=12;report['unresolved_placeholders']=0
    seeds=rows(F/'scores/seed_scores.csv');by={}
    for r in seeds:by.setdefault(r['method'],{})[r['seed']]=float(r['four_world_rectangle_MSE'])
    diff=np.array([by['OLS_high'][s]-by['OLS_low'][s] for s in sorted(by['OLS_high'])])
    mean=diff.mean();sem=diff.std(ddof=1)/math.sqrt(len(diff));ci=[mean-student_t.ppf(.975,999)*sem,mean+student_t.ppf(.975,999)*sem]
    ref=rows(F/'scores/paired_comparisons.csv')[0]
    assert max(abs(mean-float(ref['paired_risk_difference'])),abs(ci[0]-float(ref['MC95_low'])),abs(ci[1]-float(ref['MC95_high'])))<1e-12
    risk_reduction=100*(1-statistics.mean(by['OLS_high'].values())/statistics.mean(by['OLS_low'].values()))
    report['formal_recalculation']={'paired_difference':float(mean),'MC95':ci,'risk_reduction_percent':risk_reduction,'status':'PASS'}
    fold=rows(P/'Step2_ED5_Folds.csv');sse={}
    for name in ['M0','M1','R']:
        rr=[r for r in fold if r['seed']=='0' and r['model']==name]
        sse[name]={k:sum(float(r[k]) for r in rr) for k in ['C_SSE','response_SSE']}
    matched={r['model']:r for r in rows(P/'Step2_Fig5_Matched55.csv')}
    gains={}
    for name in ['M1','R']:
        cg=1-sse[name]['C_SSE']/sse['M0']['C_SSE'];rg=1-sse[name]['response_SSE']/sse['M0']['response_SSE']
        assert abs(cg-float(matched[name]['C_gain_vs_M0']))<1e-12
        assert abs(rg-float(matched[name]['response_gain_vs_M0']))<1e-12
        gains[name]={'rectangle_gain_percent':100*cg,'complete_gain_percent':100*rg}
    report['RNA_pooled_gains']=gains
    choices=rows(G/'CHOICES.csv');allg=rows(G/'SUMMARY.csv');assert len(choices)==420 and len(allg)==15
    largest=0
    for r in allg:
        cc=[c for c in choices if c['model']==r['model']];assert len(cc)==28
        assert sum(int(c['hit']) for c in cc)==int(r['hits'])
        for column,src_col in [('mean_U','actual_U'),('mean_regret','regret'),('half0_U','half0_U'),('half1_U','half1_U')]:
            error=abs(statistics.mean(float(c[src_col]) for c in cc)-float(r[column]));largest=max(largest,error);assert error<1e-12
    rg={r['group']:r['pair'] for r in choices if r['model']=='R_G'}
    pr={r['group']:r['pair'] for r in choices if r['model']=='PRESAGE_Combo_seed24'};assert rg==pr
    groups=js(G/'GROUPS.json');assert len(groups)==28 and sum(g['n']>=4 for g in groups)==5
    assert len(set(c['physical_anchor'] for c in choices))==8
    report['GSE_choice_recalculation']={'summary_rows':15,'choice_rows':420,'max_absolute_difference':largest,'RG_PRESAGE_choices_identical':True,'status':'PASS'}
    jobs=[]
    for kind in ['gears','scouter']:
        for seed in range(3):
            for f in range(5):
                d=js(G/f'runs/{kind}_seed{seed}/fold_{f}/SELECTION.json');sel=d.get('selected',d)
                run=js(G/f'runs/{kind}_seed{seed}/fold_{f}/refit/RUN.json')
                assert sel['status']==run['status']=='COMPLETE'
                assert sel['selected_epoch']==run['epochs']
                if kind=='scouter':assert not sel['extended']
                jobs.append({'model':kind,'seed':seed,'fold':f,'epochs':run['epochs'],'schedule':sel.get('schedule')})
    report['neural_jobs_checked']=len(jobs);report['neural_selection_refit_match']='PASS'
    with (ROOT/'NEURAL_TRAINING_COMPLETION.csv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=jobs[0].keys());w.writeheader();w.writerows(jobs)
    report['main_sha256']=sha(ROOT/'WitnessCell_v24_Main_Methods.md')
    report['SI_sha256']=sha(ROOT/'WitnessCell_v24_Supplementary_Information.md')
    report['SI_words_including_tables']=len(re.findall(r"\b[\w]+(?:[-–'][\w]+)*\b",text))
    report['scope']='Document/source verification and saved-result recalculation only; no model retraining.'
    (ROOT/'SI_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path);ap.add_argument('--bridge',type=Path);ap.add_argument('--main',type=Path);args=ap.parse_args()
    if args.repo:
        assert args.bridge and args.main
        snapshot(args.repo.resolve(),args.bridge.resolve(),args.main.resolve())
    text=build();verify(text)
