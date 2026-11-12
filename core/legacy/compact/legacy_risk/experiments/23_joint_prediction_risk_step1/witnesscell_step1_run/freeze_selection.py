"""Freeze all delivery identities from scores before attaching test labels."""
import csv,json,hashlib,time
from decimal import Decimal
from pathlib import Path
HERE=Path(__file__).resolve().parent
PID='WitnessCell-step1-v1-20261017'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
allfreeze=json.loads((HERE/'all_risks_frozen.json').read_text())
files={}
for name,digest in allfreeze['outputs'].items():
    d=HERE/'results'/name;p=d/'frozen_query_outputs.csv';assert sha(p)==digest
    assert not (d/'scored_queries.csv').exists(),'Already labelled: cannot claim a new pre-label freeze'
    with p.open() as f:rows=list(csv.DictReader(f))
    flags=[]
    for cov in [.3,.5,.7,.9]:
        k=int(Decimal(str(cov))*len(rows))
        for method,column in [('S','risk_S'),('B','risk_B'),('effect','effect_score'),('variance','variance_score')]:
            order=sorted(range(len(rows)),key=lambda i:(float(rows[i][column]),hashlib.sha256((PID+'|'+rows[i]['dataset']+'|'+rows[i]['pair']).encode()).hexdigest(),rows[i]['pair']))
            accepted=set(order[:k])
            for i,row in enumerate(rows):flags.append({'dataset':row['dataset'],'seed':row['seed'],'pair':row['pair'],'coverage':cov,'method':method,'accepted':int(i in accepted),'quota':k,'candidate_count':len(rows)})
    out=d/'frozen_selections.csv'
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(flags[0]));w.writeheader();w.writerows(flags)
    files[name]=sha(out)
(HERE/'all_selections_frozen.json').write_text(json.dumps({'frozen_at_unix':time.time(),'files':files,'test_labels_scored':False},indent=2)+'\n')
print(json.dumps({'selections_frozen':len(files),'test_labels_scored':False}))
