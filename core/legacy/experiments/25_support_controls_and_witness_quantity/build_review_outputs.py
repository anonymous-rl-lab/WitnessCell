#!/usr/bin/env python3
"""Descriptive exports and figures from completed, unchanged follow-up results."""
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR', '/tmp/witnesscell_followup_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
SEEN = json.loads((ROOT/'seen/evaluation/results.json').read_text())
WITNESS = json.loads((ROOT/'wessels_witness/statistics/summary.json').read_text())


def read_csv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    rows = list(rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ensure_context(src, dst):
    if src.is_file():
        shutil.copyfile(src, dst)
    elif not dst.is_file():
        raise FileNotFoundError('Neither the original nor the packaged context file exists: '+str(dst))


def seen_result(dataset, scope='all_queries', comparison='S/C', coverage=.5, loss='loss_mse'):
    return next(r for r in SEEN['results'] if (r['dataset'],r['scope'],r['comparison'],r['coverage'],r['loss'])
                ==(dataset,scope,comparison,coverage,loss))


def pct_ci(gain, ci):
    return f'{gain*100:+.2f}% [{ci[0]*100:+.2f}, {ci[1]*100:+.2f}]'


def main():
    # Historical mean-prediction evidence is copied and explicitly labelled as reused.
    source=ROOT.parent/'v19_inputs/step1_unpacked/witnesscell_step1_run'
    legacy=ROOT/'original_step1_context'; legacy.mkdir(exist_ok=True)
    inventory=[]
    for name in ('baseline_summary.json','baseline_summary.csv','baseline_seed_results.csv'):
        src=source/'evaluation'/name; dst=legacy/name
        ensure_context(src,dst)
        inventory.append({'source':str(src),'copy':str(dst.relative_to(ROOT)),'sha256':sha(dst),
                          'status':'original_Step1_result_reused_not_a_new_experiment'})
    for ds in ('Norman','Wessels'):
        for seed in (1,2,3):
            src=source/'results'/f'{ds}_seed{seed}'/'baseline_losses.csv'
            dst=legacy/f'{ds}_seed{seed}_baseline_losses.csv'
            ensure_context(src,dst)
            inventory.append({'source':str(src),'copy':str(dst.relative_to(ROOT)),'sha256':sha(dst),
                              'status':'original_Step1_per_query_losses'})
    if not (legacy/'SOURCE_INVENTORY.json').exists():
        (legacy/'SOURCE_INVENTORY.json').write_text(json.dumps(inventory,indent=2)+'\n')

    table=[]
    for ds in ('Norman','Wessels'):
        for scope in ('all_queries','seen1'):
            for comp in ('S/B','S/R','S/C','S/C0'):
                r=seen_result(ds,scope,comp)
                table.append({'dataset':ds,'scope':scope,'comparison':comp,'coverage':.5,
                    'gain_percent':100*r['gain'],'ci95_low_percent':100*r['gain_ci95'][0],
                    'ci95_high_percent':100*r['gain_ci95'][1],
                    'query_seed_units':sum(x['n'] for x in r['seed_counts']),
                    'accepted_query_seed_units':sum(x['quota'] for x in r['seed_counts']),
                    'valid_bootstrap':r['valid_bootstrap'],'invalid_bootstrap':r['invalid_bootstrap']})
    write_csv(ROOT/'tables/seen_primary_comparisons.csv',table)

    # Descriptive distributions use existing frozen decisions; there is no new selection model.
    distributions=[]; allocations=[]; overlaps=[]
    for ds in ('Norman','Wessels'):
        for seed in (1,2,3):
            folder=ROOT/'seen/results'/f'{ds}_seed{seed}'
            rows={r['pair']:r for r in read_csv(folder/'new_scored_queries.csv')}
            decisions=read_csv(folder/'new_frozen_decisions.csv')
            for scope in ('all_queries','seen1'):
                chosen={}
                for method in ('B','R','C','C0','S'):
                    rr=[r for r in decisions if r['scope']==scope and r['method']==method
                        and float(r['coverage'])==.5]
                    chosen[method]={r['pair'] for r in rr if int(r['accepted'])}
                    for accepted in (0,1):
                        group=[rows[r['pair']] for r in rr if int(r['accepted'])==accepted]
                        vals=np.array([float(r['true_effect']) for r in group])
                        q=np.quantile(vals,[0,.25,.5,.75,1])
                        distributions.append({'dataset':ds,'seed':seed,'scope':scope,'method':method,
                            'accepted':accepted,'n':len(vals),'mean_true_effect':float(vals.mean()),
                            'min':q[0],'q25':q[1],'median':q[2],'q75':q[3],'max':q[4],
                            'definition':'original_Step1_non_target_true_effect; descriptive, not independent repeats'})
                    if scope=='all_queries':
                        for seen in (0,1,2):
                            allocations.append({'dataset':ds,'seed':seed,'method':method,'seen':seen,
                                'candidate_count':sum(int(r['seen'])==seen for r in rows.values()),
                                'accepted_count':sum(int(rows[p]['seen'])==seen for p in chosen[method])})
                for method in ('B','R','C','C0'):
                    common=len(chosen['S']&chosen[method]); union=len(chosen['S']|chosen[method])
                    overlaps.append({'dataset':ds,'seed':seed,'scope':scope,'comparison':'S/'+method,
                        'quota':len(chosen['S']),'intersection':common,'union':union,'jaccard':common/union,
                        'S_only':len(chosen['S']-chosen[method]),'comparator_only':len(chosen[method]-chosen['S'])})
    write_csv(ROOT/'tables/accepted_rejected_true_effect_distributions.csv',distributions)
    write_csv(ROOT/'tables/acceptance_overlap.csv',overlaps)
    allocation_summary=[]
    for ds in ('Norman','Wessels'):
        for method in ('B','R','C','C0','S'):
            row={'dataset':ds,'method':method}
            for seen in (0,1,2):
                r=[x for x in allocations if x['dataset']==ds and x['method']==method and x['seen']==seen]
                row[f'seen{seen}_accepted']=sum(x['accepted_count'] for x in r)
                row[f'seen{seen}_candidates']=sum(x['candidate_count'] for x in r)
            allocation_summary.append(row)
    write_csv(ROOT/'tables/seen_allocation_pooled_descriptive.csv',allocation_summary)

    figdir=ROOT/'figures'; figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12.8,7.0),gridspec_kw={'width_ratios':[1.05,1]})
    fig.subplots_adjust(left=.16,right=.98,bottom=.34,top=.81,wspace=.32)
    colors={'Norman':'#A65532','Wessels':'#176B86'}
    locations=[3,2,1,0]
    configurations=[('Norman','all_queries'),('Norman','seen1'),('Wessels','all_queries'),('Wessels','seen1')]
    labels=[]
    for y,(ds,scope) in zip(locations,configurations):
        r=seen_result(ds,scope); val=r['gain']*100; lo,hi=np.array(r['gain_ci95'])*100
        ax.errorbar(val,y,xerr=[[val-lo],[hi-val]],fmt='s' if scope=='seen1' else 'o',
                    color=colors[ds],markersize=7,capsize=4,elinewidth=1.8)
        ax.text(val,y+.19,f'{val:+.1f}%',ha='center',color=colors[ds],fontsize=10)
        labels.append(ds+'\n'+('all queries' if scope=='all_queries' else 'within seen1'))
    ax.axvline(0,color='#6B7280',linewidth=1,linestyle='--')
    ax.set(yticks=locations,yticklabels=labels,xlim=(-23,31),ylim=(-.5,3.65),
           xlabel='S versus seen-category control C\nNon-target MSE reduction (%)')
    ax.set_title('a  Increment beyond seen categories',loc='left',fontweight='bold',pad=16)
    ax.grid(axis='x',alpha=.18)
    for metric,label,color,marker in [('loss_mse','MSE (primary)','#176B86','o'),
                                      ('loss_wmse','WMSE (secondary)','#9B6B22','s')]:
        rows=[r for r in WITNESS['fraction_results'] if r['metric']==metric and r['reference']=='factorized_prediction']
        x=np.array([r['fraction'] for r in rows])*100
        gain=np.array([r['gain'] for r in rows])*100
        ci=np.array([r['gain_ci95'] for r in rows])*100
        bx.plot(x,gain,color=color,marker=marker,label=label,linewidth=2)
        bx.fill_between(x,ci[:,0],ci[:,1],color=color,alpha=.12)
        if metric=='loss_mse':
            for a,b in zip(x,gain): bx.annotate(f'{b:.1f}%',(a,b),xytext=(0,-17),textcoords='offset points',ha='center',color=color)
    bx.set(xlim=(18,106),ylim=(10,52),xticks=[25,50,100],xlabel='Available residual witness library (%)',
           ylabel='Non-target error reduction versus\nfixed factorized backbone (%)')
    bx.set_title('b  Wessels witness-library intervention',loc='left',fontweight='bold',pad=16)
    bx.grid(axis='y',alpha=.18); bx.legend(frameon=False,loc='upper left',fontsize=9)
    fig.suptitle('WitnessCell: support adds conditional value',x=.55,y=.96,fontsize=17,fontweight='bold')
    fig.text(.55,.885,'Frozen final predictor for selection; frozen endpoints and parameters for the witness intervention',ha='center',fontsize=10)
    fig.text(.16,.07,'Points / bands: 95% pair-cluster bootstrap intervals (5,000 draws). Retrospective analyses.\n'
             'Norman S/C primary comparison is inconclusive. Wessels effect-quartile control also spans zero.\n'
             'Witness intervention: 100% minus 25% MSE gain = 14.15 percentage points [12.42, 16.06].',fontsize=9,linespacing=1.6)
    for ext in ('png','pdf','svg'): fig.savefig(figdir/f'step1_followup_summary.{ext}',dpi=180)
    plt.close(fig)

    core={'seen':table,'effect_control':SEEN['effect_quartile_control'],
          'witness_fraction_results':WITNESS['fraction_results'],'witness_contrasts':WITNESS['contrasts'],
          'historical_mean_prediction_results':[r for r in read_csv(legacy/'baseline_summary.csv')
             if r['pool']=='all_queries' and r['comparator']=='shared_endpoint_additive'],
          'primary_conclusion':'Norman S/C did not establish incremental benefit beyond seen; Wessels supports a bounded conditional increment. Fixed-state Wessels witness-quantity effect is positive.',
          'source_files':{'seen':'seen/evaluation/results.json','witness':'wessels_witness/statistics/summary.json',
                          'baseline':'original_step1_context/baseline_summary.csv'}}
    (ROOT/'KEY_RESULTS.json').write_text(json.dumps(core,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'DESCRIPTIVE_EXPORTS_COMPLETE','seen_comparisons':len(table),
                      'true_effect_distribution_rows':len(distributions),'figure':'figures/step1_followup_summary.png'}))


if __name__=='__main__': main()
