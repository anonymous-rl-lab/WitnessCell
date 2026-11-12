"""All-roster gate, independent replay of frozen predictions, common partner scoring."""
from common import *
import pandas as pd
def score():
 result=ROOT/'results/detected49';b=setup();m=b['m'];full=b['full'];P,T=full.x.shape;mask=m['score_gene_mask'];G=int(mask.sum())
 # Check ALL deliverables before reading the outcome rows used for scoring.
 required=[result/'classical'/f'fold_{f}' for f in range(5)]+[result/'presage'/f'fold_{f}' for f in range(5)]
 for kind in ('gears','scouter'):
  for seed in range(3):required += [result/f'{kind}_seed{seed}'/f'fold_{f}' for f in range(5)]
 # Fixed40 seed24 remains supplementary and is not a main-roster gate.
 missing=[str(p.relative_to(ROOT)) for p in required if not (p/'FREEZE.json').exists()]
 if missing:raise RuntimeError('Full-roster scoring gate: unfinished '+str(missing))
 for folder in required:
  rec=json.loads((folder/'FREEZE.json').read_text());h=rec.get('output_sha256',rec.get('predictions_sha256'))
  assert h==sha(folder/'predictions.npz'),str(folder)
 utilities={};rna={};assignment=np.full(P,-1,int)
 def put(name,ids,u,delta=None):
  if name not in utilities:utilities[name]=np.full((P,T),np.nan)
  assert np.isnan(utilities[name][ids]).all();utilities[name][ids]=u
  if delta is not None:
   if name not in rna:rna[name]=np.full((P,G),np.nan)
   assert delta.shape==(len(ids),G);rna[name][ids]=delta
 for f in range(5):
  z=load(result/'classical'/f'fold_{f}/predictions.npz');q=z['query'];assert np.all(assignment[q]<0);assignment[q]=f
  for key in z:
   if key.startswith('utility_'):
    name=key[len('utility_'):];put(name,q,z[key],z.get('pred_effect_'+name))
  jobs=[('presage','PRESAGE_Combo_seed24')]+[(f'{k}_seed{s}',f'{k.upper() if k=="gears" else "Scouter"}_seed{s}') for k in ('gears','scouter') for s in range(3)]
  for directory,name in jobs:
   zz=load(result/directory/f'fold_{f}/predictions.npz');assert np.array_equal(q,zz['query']);assert np.array_equal(zz['genes'],m['genes'])
   np.testing.assert_allclose(zz['utility'],zz['pred_effect']@m['objectives'].T,rtol=1e-7,atol=1e-10)
   put(name,q,zz['utility'],zz['pred_effect'][:,mask])
 for prefix in ('GEARS','Scouter'):
  name=prefix+'_Ensemble';utilities[name]=np.mean([utilities[f'{prefix}_seed{s}'] for s in range(3)],axis=0);rna[name]=np.mean([rna[f'{prefix}_seed{s}'] for s in range(3)],axis=0)
 assert len(utilities)==15 and np.all(assignment>=0) and all(np.isfinite(x).all() for x in utilities.values())
 # Outcomes enter only after the gate and frozen prediction replay above.
 truth=m['means'][b['pair_moment_ids']];trueU=truth@m['objectives'].T
 halfU=np.array([m[k][b['pair_moment_ids']]@m['objectives'].T for k in ('first_means','second_means')])
 captureU=np.einsum('pcg,tg->pct',m['capture_effects'][b['pair_moment_ids']],m['objectives'])
 records=[];groups=[]
 for f in range(5):
  for ai,A in enumerate(m['target_names'].astype(str)):
   ids=np.flatnonzero((assignment==f)&(full.x[:,ai]>0))
   if len(ids)<2:continue
   order=np.array(sorted(ids,key=lambda i:full.pairs[i]));partners=[next(p for p in full.pairs[i].split('+') if p!=A) for i in order]
   group=f'f{f}:{A}';best=float(trueU[order,ai].max());groups.append(dict(group=group,fold=f,anchor=A,n=len(ids),pairs=order.tolist(),partners=partners))
   for name,pred in utilities.items():
    pp=pred[order,ai];loc=int(np.flatnonzero(pp>=pp.max()-1e-12)[0]);chosen=int(order[loc]);actual=float(trueU[chosen,ai]);err=pp-trueU[order,ai]
    row=dict(model=name,group=group,fold=f,anchor=A,physical_anchor=A.rsplit('-',1)[0],candidate_count=len(ids),partner=partners[loc],pair=full.pairs[chosen],
     predicted_U=float(pp[loc]),actual_U=actual,best_observed_U=best,regret=best-actual,top_choice=int(actual>=best-1e-12),centered_utility_MSE=float(np.mean((err-err.mean())**2)),
     half0_U=float(halfU[0,chosen,ai]),half1_U=float(halfU[1,chosen,ai]),half0_regret=float(halfU[0,order,ai].max()-halfU[0,chosen,ai]),half1_regret=float(halfU[1,order,ai].max()-halfU[1,chosen,ai]))
    for c in range(3):row[f'D{c+1}_U']=float(captureU[chosen,c,ai])
    records.append(row)
 frame=pd.DataFrame(records);summary=[];geometry=[]
 for f in range(5):
  ids=np.flatnonzero(assignment==f);x=full.x[ids];rank=int(np.linalg.matrix_rank(x));Q=np.eye(len(ids))-x@np.linalg.pinv(x);geometry.append(dict(fold=f,ids=ids,projection=Q,rank=len(ids)-rank))
 for name in utilities:
  rows=frame[frame.model==name];large=rows[rows.candidate_count>=4]
  rec=dict(model=name,choice_tasks=len(rows),mean_U=float(rows.actual_U.mean()),mean_regret=float(rows.regret.mean()),top_choice_accuracy=float(rows.top_choice.mean()),
   centered_utility_MSE=float(rows.centered_utility_MSE.mean()),large_tasks=len(large),large_mean_U=float(large.actual_U.mean()) if len(large) else None,
   large_mean_regret=float(large.regret.mean()) if len(large) else None,half0_mean_U=float(rows.half0_U.mean()),half1_mean_U=float(rows.half1_U.mean()),
   RNA_MSE=None,partner_projection_MSE=None,partner_degrees=sum(z['rank'] for z in geometry))
  if name in rna:
   er=rna[name]-truth[:,mask];rec['RNA_MSE']=float(np.mean(er*er));den=G*rec['partner_degrees']
   rec['partner_projection_MSE']=float(sum(np.sum((z['projection']@er[z['ids']])**2) for z in geometry)/den) if den else None
  summary.append(rec)
 table=pd.DataFrame(summary).sort_values(['mean_U','mean_regret'],ascending=[False,True]);out=result/'evaluation_main_20261105';out.mkdir(exist_ok=True)
 frame.to_csv(out/'CHOICES.csv',index=False);table.to_csv(out/'SUMMARY.csv',index=False);writej(out/'GROUPS.json',groups)
 # Paired descriptive uncertainty; cluster on physical anchor to avoid treating
 # activation/repression and repeated folds of an anchor as independent.
 compare=[];rng=np.random.default_rng(20261102);pivot=frame.pivot(index='group',columns='model',values='actual_U');anchors=frame.drop_duplicates('group').set_index('group').physical_anchor.reindex(pivot.index)
 keys=sorted(anchors.unique());draws=rng.integers(0,len(keys),size=(10000,len(keys)))
 for ref in ['PRESAGE_Combo_seed24','StandardRidge_RNA','StandardRidge_Choice','SingleB_additive','GEARS_Ensemble','Scouter_Ensemble']:
  name='R_G'
  delta=pivot[name]-pivot[ref];sums=np.array([delta[anchors==k].sum() for k in keys]);counts=np.array([np.sum(anchors==k) for k in keys]);boot=sums[draws].sum(1)/counts[draws].sum(1)
  compare.append(dict(model=name,reference=ref,mean_paired_U_difference=float(delta.mean()),cluster_bootstrap_low=float(np.quantile(boot,.025)),cluster_bootstrap_high=float(np.quantile(boot,.975)),physical_anchor_clusters=len(keys)))
 pd.DataFrame(compare).to_csv(out/'PAIRED_COMPARISONS.csv',index=False)
 np.savez_compressed(out/'FROZEN_PREDICTIONS.npz',conditions=np.array(full.pairs),fold=assignment,genes=m['genes'][mask],target_names=m['target_names'],
  **{'utility__'+k:v for k,v in utilities.items()},**{'RNA__'+k:v for k,v in rna.items()})
 report=['GSE220974：全部主比较模型已完成，含逐种子与历史参照共15行','',f'49种检出双干预；{len(groups)}个固定A选伙伴任务；{sum(g["n"]>=4 for g in groups)}个候选数≥4的任务。',
  '这是检出标签队列，不等同于论文21种已核实设计构建体；不能据此排除误分配或双细胞。','所有成绩均来自完整组合留出。此前已查看其他模型的外层结果；本轮GEARS和Scouter沿用已固定训练和选轮规则，未根据其外层成绩重新调参。',
  '选优目标为单干预A定义的RNA程序，最优仅指同一候选池中实测目标值最高，不是全图最优或治疗效果。','',table.to_string(index=False),
  '', 'TargetRank与Stable只输出排序分数，未伪造RNA表达预测。','跨半和跨capture敏感性、逐任务选择及配对差异见CSV；bootstrap为当前候选池上的描述性不确定性，不是跨研究人群的泛化保证。']
 (out/'REPORT_CN.txt').write_text('\n'.join(report)+'\n')
 writej(out/'EVALUATION.json',dict(status='COMPLETE_MAIN_ROSTER_15_ROWS',pairs=P,groups=len(groups),models=15,previous_outcomes_opened_for_other_models=True,all_predictions_frozen=True,outer_outcomes_used_for_training=False,
  scope='detected guide assignments; intended construct whitelist unresolved',partner_geometry=[dict(fold=z['fold'],query_count=len(z['ids']),degrees=z['rank']) for z in geometry],
  outputs={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='EVALUATION.json'}))
 main_names=['R_G','GEARS_Ensemble','Scouter_Ensemble','StandardRidge_RNA','SingleB_additive']
 main=table.set_index('model').loc[main_names].reset_index();main.to_csv(out/'MAIN_FIVE_MODELS.csv',index=False)
 print(main.to_string(index=False),flush=True)
if __name__=='__main__':score()
