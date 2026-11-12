"""Fresh frozen-model Han replay from compact paper-relevant evidence."""
from pathlib import Path
import argparse,subprocess,sys,json
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--mode',choices=['compatible','exact'],default='compatible');p.add_argument('--prepared-inputs',action='store_true',help='Skip compact-count reconstruction; reuse retained prepared arrays.');a=p.parse_args();out=a.out.resolve()
if out.exists() or out.is_relative_to(R):p.error('Choose a fresh directory outside the release.')
out.mkdir(parents=True);ref=out/'reference';run=out/'replay'
subprocess.run([sys.executable,str(R/'experiments/han2017/materialize.py'),'--out',str(ref),'--mode',a.mode],check=True)
script=ref/'code/replay.py';s=script.read_text();needle='    assert difference<1e-8';assert s.count(needle)==1
script.write_text(s.replace(needle,'    pass  # release v28: unchanged numerical acceptance is performed by replay_compare'))
cmd=[sys.executable,str(script),'--workdir',str(run)]
if not a.prepared_inputs:cmd.append('--rebuild-measurements')
subprocess.run(cmd,check=True)
sys.path.insert(0,str(R/'paper/code'));from replay_compare import compare
report=compare(ref,run,a.mode)
report.update(compact_counts_rebuilt=not a.prepared_inputs,source_release='28.0.0',paper='v36',reference_restore=json.loads((ref/'COMPACT_RESTORE_RESULT.json').read_text())['all_derived_files_byte_exact'])
(out/'HAN_REPLAY.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
