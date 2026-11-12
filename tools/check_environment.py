"""Report the replay environment; never install or modify dependencies."""
from pathlib import Path
import argparse,importlib.metadata as m,json,platform,sys
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--strict',action='store_true');ap.add_argument('--out',type=Path);a=ap.parse_args()
rows=[]
for line in (R/'requirements-tested.txt').read_text().splitlines():
 line=line.strip()
 if not line or line.startswith('#'):continue
 name,expected=line.split('==')
 try:version=m.version(name)
 except m.PackageNotFoundError:version=None
 rows.append(dict(package=name,expected=expected,installed=version,exact=version==expected))
report=dict(python=sys.version,platform=platform.platform(),packages=rows,all_present=all(x['installed'] for x in rows),pinned_versions_match=all(x['exact'] for x in rows),python_312=sys.version_info[:2]==(3,12))
if a.out:
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(0 if report['all_present'] and (not a.strict or report['pinned_versions_match'] and report['python_312']) else 1)
