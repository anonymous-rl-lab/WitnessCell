#!/usr/bin/env python3
"""Apply literal official GEARS pandas top-(k+1) graph preprocessing."""
from pathlib import Path
import hashlib,json,time
import pandas as pd
P=Path(__file__).resolve().parent;d=P/'go_essential_all';src=d/'go_essential_all.csv';dst=d/'go_essential_all_top21.csv'
t0=time.monotonic()
f=pd.read_csv(src,dtype={'source':'category','target':'category','importance':'float64'})
print('read',len(f),'rows',time.monotonic()-t0,flush=True)
out=f.groupby('target',observed=True).apply(lambda x:x.nlargest(21,['importance'])).reset_index(drop=True)
# Categories only change storage, not group ordering or each group's source-row order.
out['source']=out['source'].astype(str);out['target']=out['target'].astype(str)
out.to_csv(dst.with_suffix('.tmp'),index=False);dst.with_suffix('.tmp').replace(dst)
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
receipt={'source_sha256':h(src),'source_rows':len(f),'source_bytes':src.stat().st_size,'rule':'Literal official pandas groupby(target).apply(lambda x:x.nlargest(21,[importance])).reset_index(drop=True)','target_count':int(f.target.nunique()),'compact_rows':len(out),'compact_sha256':h(dst),'compact_bytes':dst.stat().st_size,'exact_implementation':'Official expression applied directly to complete reconstructed graph. Categorical node strings reduce RAM, preserving group and row ordering.','wall_seconds':time.monotonic()-t0}
(P/'go_graph_top21_provenance.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
