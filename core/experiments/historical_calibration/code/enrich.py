import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,re
R=Path(__file__).resolve().parents[1];H=R.parent/'gain_regression_history';S=R.parent/'gain_regression_sources';OLD=R.parent/'partner_gain_regression_20261108';O=R/'results'
d=pd.read_csv(R/'inputs/HISTORICAL_INPUT.csv');d['record_id']=np.arange(len(d));ROWS=[];CHECKS=[];EVIDENCE={};INNER_GROWTH={}
def js(p):return json.loads(Path(p).read_text())
def xof(e,n):
 x=np.zeros((len(e),n));x[np.arange(len(e))[:,None],e]=1;return x
def proj(e,n,y):
 x=xof(e,n);return y-x@np.linalg.lstsq(x,y,rcond=None)[0]
def edges(names):
 pairs=[str(s).split('+') for s in names];nodes=sorted(set(sum(pairs,[])));ix={s:i for i,s in enumerate(nodes)}
 return np.array([[ix[a],ix[b]] for a,b in pairs]),len(nodes)
def summary(y,base,pred,e,n):
 pe=proj(e,n,y-base);pc=proj(e,n,pred-base);be=y-base;c=pred-base;b=np.sum(pe*pe);full=np.sum(be*be)
 return dict(p_base=float(b),p_loss=float(np.sum((pe-pc)**2)),f_base=float(full),f_loss=float(np.sum((be-c)**2)),reward=float(2*np.sum(pe*pc)),cost=float(np.sum(pc*pc)),inner_dims=len(e)-np.linalg.matrix_rank(xof(e,n)))
def combine(stats):
 sums={k:sum(s[k] for s in stats) for k in stats[0]};return dict(inner_partner_gain=100*(1-sums['p_loss']/sums['p_base']),inner_full_gain=100*(1-sums['f_loss']/sums['f_base']),inner_reward=100*sums['reward']/sums['p_base'],inner_cost=100*sums['cost']/sums['p_base'],inner_dims=sums['inner_dims'],inner_fold_gain_sd=float(np.std([100*(1-z['p_loss']/z['p_base']) for z in stats])),inner_folds=len(stats),inner_metric='foldwise_endpoint_projection',validation_reference='M0')
def put(row,extra,origin,pool):
 rec={**row,'zero_to_reference':False,**extra,'inner_origin':str(origin),'pool_id':pool};ROWS.append(rec)
 EVIDENCE[str(origin)]=hashlib.sha256(Path(origin).read_bytes()).hexdigest()
def replay(models,singles,e):
 scale=np.std(singles,axis=0);scale[scale==0]=1;i,j=e.T;out=[]
 for t,m in enumerate(models):
  a,ep,db,dc=[np.asarray(m[k]) for k in ['a','e','d_beta','d_c']]
  v=.5*(a[i]+a[j]+(m['beta']+db[i])*a[j]+(m['beta']+db[j])*a[i]+(m['c']+dc[i])*a[j]**2+(m['c']+dc[j])*a[i]**2)+m['e0']+ep[i]+ep[j];out.append(v*scale[t])
 return np.column_stack(out)
# Recover functional inner rank-one predictions from stored u/lambda; verify full SSE.
feat=d[d.variant=='raw'].drop_duplicates('feature_id');pilot={v['label']:v for v in js(OLD/'inputs/metadata.json')}
for _,fr in feat.iterrows():
 key=fr.feature_id
 if key in pilot:
  meta=pilot[key];a=np.load(OLD/'inputs'/meta['train_file']);e=a['edges'];y=a['y'];n=len(a['genes']);folds=a['fold'];oof=a['M0_oof'];cell=meta['cell']
  if cell=='K562':md=S/'horlbeck/WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/K562/cases'/f"outer{meta['outer_fold']}"/'results'/meta['budget']
  elif cell=='Jurkat':md=S/'horlbeck/WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/Jurkat/results/high'
  else:md=S/'WitnessCell_Shen2017_Formal_20261108/WitnessCell_Shen2017_Formal_20261108/runs/formal'/cell/'results/high'
 elif key.startswith('K562_layout_'):
  mt=re.fullmatch(r'K562_layout_f(\d)_s(\d)_(.*)',key);fold,seed,arm=mt.groups();b=H/'bridge_restored/WitnessCell_NM_NCS_Bridge_20261031/growth_run';md=b/'training_v2'/f'fold_{fold}_seed_{seed}'/'results'/arm
  a=np.load(S/'horlbeck/WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/K562/cases'/f'outer{fold}'/'inputs/training.npz');la=np.load(b/'layouts_v2'/f'fold_{fold}_seed_{seed}'/f'{arm}.npz');e=la['edges'];y=a['y'][la['positions']];n=len(a['genes']);folds=la['inner_fold'];oof=np.zeros_like(y)
  for ff in range(3):
   ix=folds==ff;oof[ix]=replay(js(md/f'inner{ff}_baseline.json'),a['singles'],e[ix])
 else:continue
 v=js(md/'inner_summary.json');sel=js(md/'selection.json');keymodel=sel['rank1']['best_regularized'];z=next(t for t in v['grid'] if t['key']==keymodel);pred=oof.copy();ranks=[]
 for diag in z['diagnostics']:
  f=int(diag['inner_fold']);ix=folds==f;u=np.asarray(diag['u']);c=float(diag['lambda'])*u[e[ix,0]]*u[e[ix,1]];support=np.bincount(e[~ix].ravel(),minlength=n)>0;c[~(support[e[ix,0]]&support[e[ix,1]])]=0;pred[ix]+=c[:,None];ranks.append(diag['numerical_rank'])
 delta=float(np.max(abs(np.sum((y-pred)**2,axis=0)-np.asarray(z['sse']))));assert delta<1e-8,(key,delta);CHECKS.append(dict(key=key,inner_full_sse_replay_error=delta))
 extra=combine([summary(y[folds==f],oof[folds==f],pred[folds==f],e[folds==f],n) for f in range(3)])
 extra.update(amplitude=1.,log_regularization=np.log10(z['nu']),model_rank=float(np.mean(ranks)),effective_df=np.nan,inner_status='archived_regularized',chosen_by_full=int(sel['rank1']['selected']!='G0'))
 # Compact inner observations permit independent calibration replay without source archives.
 for suffix,value in {'y':y,'M0':oof,'regularized':pred,'edges':e,'folds':folds,'nodes':np.array(n),'selected_zero':np.array(sel['rank1']['selected']=='G0')}.items():INNER_GROWTH[key+'__'+suffix]=value
 # Retain selected outputs and regularized diagnostic outputs as distinct candidate rows.
 candidates=d[(d.feature_id==key)&(d.variant=='raw')&d.model.isin(['rank1','rank1_selected','rank1_regularized'])].copy()
 if key in pilot:
  aout=np.load(md/'predictions.npz');truth=np.load(OLD/'inputs'/pilot[key]['outcome_file']);sse=float(np.sum(proj(aout['edges'],n,truth['y']-aout['rank1_regularized'])**2));den=float(np.sum(proj(aout['edges'],n,truth['y']-aout['M0'])**2));new=fr.to_dict();new.update(record_id=-len(CHECKS),model='rank1_regularized',partner_sse=sse,base_partner_sse=den,gain_pp=100*(1-sse/den),source=str(md/'predictions.npz'));candidates=pd.concat([candidates,pd.DataFrame([new])],ignore_index=True)
 for row in candidates.to_dict('records'):
  iszero=(row['model']!='rank1_regularized' and sel['rank1']['selected']=='G0');ex=dict(extra)
  if iszero:ex.update(inner_partner_gain=0.,inner_full_gain=0.,inner_reward=0.,inner_cost=0.,amplitude=0.,model_rank=0.,inner_status='zero_selected',zero_to_reference=True)
  put(row,ex,md/'selection.json','growth_'+key)
np.savez_compressed(R/'inputs/INNER_GROWTH_REPLAY.npz',**INNER_GROWTH)
print('functional matched',len(ROWS),flush=True)
# RNA R archives: saved inner prediction arrays, reconstructed endpoint-projection scores.
for path,group in d[(d.variant=='raw')&d.source.str.contains('/predictions.npz')&d.panel.eq('CaRPool')].groupby('source'):
 p=Path(path);md=p.parent;ps=sorted(md.glob('inner_selected_*.npz'))
 if len(ps)!=3 or not (md/'selection.json').exists():continue
 sel=js(md/'selection.json');by={};rstats=[]
 for pp in ps:
  a=np.load(pp);e,n=edges(a['conditions']);base=a['base'];y=a['y'];rstats.append(summary(y,base,a['R'],e,n))
  if 'model_names' in a.files:
   values=a['selected_prediction_values'];models={str(k):values[j] for k,j in zip(a['model_names'],a['model_prediction_index'])}
  else:models={k[5:]:a[k] for k in a.files if k.startswith('pred_')}
  models.update(R_original=a['R'],R=a['R'])
  for name,pred in models.items():by.setdefault(name,[]).append(summary(y,base,pred,e,n))
 for row in group.to_dict('records'):
  name=row['model']
  if name not in by or len(by[name])!=3:continue
  ex=combine(by[name]);choice=sel.get('selected',{}).get(name,{})
  ex.update(amplitude=choice.get('beta',choice.get('alpha',1.)),log_regularization=np.log10(choice.get('lam',1.) or 1.),effective_df=np.nan,model_rank=np.nan,inner_status='selected_on_same_inner_results',chosen_by_full=np.nan)
  put(row,ex,md/'selection.json',str(md))
print('R matched',len(ROWS),flush=True)
# Early RNA feature/kernel probes: exact candidate-specific inner P losses and complexity.
B=H/'WitnessCell_Dixit_Representation_Debug_20261020/WitnessCell_Dixit_Representation_Debug_20261020'
for row in d[(d.variant=='raw')&d.model.str.startswith(('core_features_','input_probe_'))].to_dict('records'):
 fold=int(row['feature_id'].rsplit('f',1)[1]);task='Wessels55' if row['panel']=='CaRPool' else 'Dixit';name=row['model'];amp=1.;lam=1.;k=None;edf=np.nan;rank=np.nan;align=np.nan
 if name.startswith('core_features_'):
  method=name[len('core_features_'):];p=B/'core_features/run1'/task/f'FOLD_{fold}.json';v=js(p);zeros=v['core_losses'];loss=v['losses'];dims=v['inner_contrast_dimensions']
  if method=='legacy_M0':continue
  if '_lambda' in method:k,la=method.split('_lambda');lam=float(la)
  else:
   choice=v['choices'][method];k=choice.get('family',method.rsplit('_',1)[0]);lam=choice.get('lambda') or 1.
  if k=='zero':amp=0.;pg=fg=0.
  else:
   pg=100*(1-loss['P'][k][str(float(lam))]/zeros['P']);fg=100*(1-loss['full'][k][str(float(lam))]/zeros['full']);au=v['outer_diagnostics']['kernels'][k];rank=au.get('effective_rank',np.nan);align=au.get('train_half_alignment',np.nan)
  scope='core_features'
 else:
  method=name[len('input_probe_'):];p=B/'input_probe/run1'/task/f'fold_{fold}/SELECTION.json';v=js(p);zeros=v['inner_zero_SSE'];loss=v['inner_SSE'];dims=v['inner_contrast_dimensions']
  if method=='legacy_M0':continue
  if method.startswith('fixed_'):k,la=method[len('fixed_'):].rsplit('_',1);lam=float(la)
  elif method in ['input_selected_P','input_selected_full','input_selected_P_gated']:
   objective='full' if method=='input_selected_full' else 'P';choice=v['joint_choices'][objective];k=choice['kernel'];lam=choice['lambda'] or 1.
   if method.endswith('_gated') and not v['gated_partner_enabled']:k='zero'
  else:k,objective=method[len('input_'):].rsplit('_',1);lam=v['per_kernel_lambda'][objective][k]
  if k=='zero':amp=0.;pg=fg=0.
  else:
   pg=100*(1-loss['P'][k][str(float(lam))]/zeros['P']);fg=100*(1-loss['full'][k][str(float(lam))]/zeros['full']);au=v['outer_kernel_audit'][k];rank=au.get('kernel_effective_rank',np.nan);edf=au.get('residual_degrees_at_lambda',{}).get(str(float(lam)),np.nan);align=au.get('training_half_covariance_alignment',np.nan)
  scope='input_probe'
 if amp==0:rank=edf=0.
 ex=dict(inner_partner_gain=pg,inner_full_gain=fg,inner_dims=dims,inner_metric='foldwise_endpoint_projection',validation_reference=row['baseline'],amplitude=amp,log_regularization=np.log10(lam),effective_df=edf,model_rank=rank,kernel_half_alignment=align,inner_status='fixed_candidate' if ('lambda' in method or method.startswith('fixed_')) else 'selected_on_same_inner_results',inner_folds=3,zero_to_reference=(amp==0))
 put(row,ex,p,scope+'_'+row['feature_id'])
print('kernel probes matched',len(ROWS),flush=True)
# Selected 12-parameter Dixit branch: same M0 comparison and inner partner score.
rr=H/'WitnessCell_Dixit_Real_Training_Test_20261026'
for row in d[(d.variant=='raw')&(d.model=='shared12_selected')].to_dict('records'):
 f=int(row['feature_id'].rsplit('f',1)[1]);p=rr/f'FOLD_{f}.json';v=js(p)
 ex=dict(inner_partner_gain=100*(1-v['inner_selected_partner_SSE']/v['inner_M0_partner_SSE']),inner_full_gain=100*(1-v['inner_selected_full_SSE']/v['inner_M0_full_SSE']),amplitude=v['selected_amplitude'],log_regularization=np.log10(v['selected_rho']),effective_df=np.nan,model_rank=np.nan,trainable_parameters=12,inner_dims=np.nan,inner_metric='foldwise_endpoint_projection',validation_reference='M0',inner_status='selected_on_same_inner_results',zero_to_reference=(v['selected_amplitude']==0))
 put(row,ex,p,'shared12_'+row['feature_id'])
# Joung v4: model-specific held-combination validation; selected and terminal both retained.
root=H/'WitnessCell_Neural_Partner_v4_TrainTest_20261025/WitnessCell_Neural_Partner_v4_TrainTest_20261025'
meta=js(root/'results/FIT_METADATA.json')
for row in d[(d.variant=='raw')&d.model.str.startswith('network_v4_')&d.panel.eq('Joung')].to_dict('records'):
 name=row['model'][len('network_v4_'):];arm,which=name.rsplit('_',2)[:2]
 # explicit suffix parsing avoids underscores inside arm names
 which='selected' if name.endswith('_selected_delta') else 'terminal';arm=name[:-len('_'+which+'_delta')]
 p=root/('energy_fix/SELECTION.json' if arm=='extra_energy' else 'results/SELECTION.json');v=js(p)
 z=v['chosen'] if arm=='extra_energy' and which=='selected' else v['per_arm'][arm] if which=='selected' else max([c for c in v['candidates'] if c['arm']==arm],key=lambda c:c['step'])
 pars=[m['trainable_parameters'] for key,m in meta.items() if key.endswith('/'+arm) and 'trainable_parameters' in m]
 ex=dict(inner_partner_gain=100*(1-z['C']),inner_full_gain=100*(1-z['F']),inner_dims=3,inner_folds=5,inner_metric='equal_fold_projection_ratios',validation_reference='M0',amplitude=float(z['step']>0),training_steps=z['step'],trainable_parameters=float(np.mean(pars)) if pars else np.nan,effective_df=np.nan,inner_status=which,zero_to_reference=(z['step']==0))
 put(row,ex,p,'Joung_v4')
print('neural matched',len(ROWS),flush=True)
# Main enriched records. Eligibility independent of outcome value.
out=pd.DataFrame(ROWS);out['inner_gain_per_dim']=out.inner_partner_gain/(1+out.inner_dims.fillna(0));out['log_inner_dims']=np.log1p(out.inner_dims);out['log_effective_df']=np.log1p(out.effective_df);out['cost_ratio']=out.inner_cost/100
out['candidate_key']=[hashlib.sha256(f"{r.feature_id}|{r.baseline}|{r.model}|{r.inner_partner_gain:.9g}|{r.inner_full_gain:.9g}|{r.amplitude:.9g}|{r.source}".encode()).hexdigest()[:20] for r in out.itertuples()]
# Use both score and inner features to identify equivalent archive aliases; no outcome is used as a regression predictor.
out['score_key']=[hashlib.sha256(f"{r.feature_id}|{r.baseline}|{r.gain_pp:.10g}|{r.inner_partner_gain:.10g}|{r.inner_full_gain:.10g}|{r.amplitude:.10g}".encode()).hexdigest()[:20] for r in out.itertuples()]
out.to_csv(O/'ENRICHED_INPUT.csv',index=False);coverage=d[['record_id','feature_id','panel','model','variant','source']].copy();coverage['inner_matched']=coverage.record_id.isin(out.record_id)
coverage['reason']=np.where(coverage.inner_matched,'matched',np.where(coverage.variant!='raw','separate_score_variant','no_matching_inner_partner_record_recovered'))
coverage.to_csv(O/'COVERAGE.csv',index=False)
(O/'INNER_REPLAY_CHECKS.json').write_text(json.dumps(CHECKS,indent=2));(O/'INNER_SOURCE_HASHES.json').write_text(json.dumps(EVIDENCE,indent=2))
print(out.groupby('panel').agg(rows=('model','size'),configs=('feature_id','nunique')).to_string(),flush=True)
print('total',len(out),'feature conditions',out.feature_id.nunique(),'source panels',out.panel.nunique(),flush=True)
