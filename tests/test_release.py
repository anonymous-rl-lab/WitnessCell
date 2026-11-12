"""Non-destructive regression tests for the release entry and coverage contract."""
from pathlib import Path
import csv,json,subprocess,sys,tempfile,unittest
R=Path(__file__).resolve().parents[1]
class ReleaseTests(unittest.TestCase):
 def cli(self,*args):
  return subprocess.run([sys.executable,'-B',str(R/'reproduce.py'),*map(str,args)],capture_output=True,text=True)
 def test_companion_mapping_has_recoverable_source(self):
  aliases={x['path']:x['canonical'] for x in json.loads((R/'provenance/ALIASES.json').read_text())}
  with (R/'provenance/V41_INTEGRATION_MAP.csv').open() as f:rows=list(csv.DictReader(f))
  self.assertEqual(len(rows),180)
  for row in rows:
   path=row['integrated_path'];source=R/aliases.get(path,path)
   self.assertTrue(source.is_file(),path)
 def test_task_inventory(self):
  r=self.cli('--list');self.assertEqual(r.returncode,0,r.stderr)
  self.assertIn('reference_sensitivity',r.stdout);self.assertIn('saturation_refit',r.stdout)
 def test_unknown_task_preserves_output(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'new';r=self.cli('--out',out,'--tasks','does_not_exist')
   self.assertEqual(r.returncode,2);self.assertFalse(out.exists())
 def test_duplicate_task_rejected(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'new';r=self.cli('--out',out,'--tasks','coverage','coverage')
   self.assertEqual(r.returncode,2);self.assertFalse(out.exists())
 def test_nested_output_rejected(self):
  r=self.cli('--out',R/'unsafe-output','--tasks','coverage');self.assertEqual(r.returncode,2);self.assertFalse((R/'unsafe-output').exists())
 def test_ancestor_output_rejected(self):
  r=self.cli('--out',R.parent,'--resume','--tasks','coverage');self.assertEqual(r.returncode,2)
 def test_existing_output_preserved(self):
  with tempfile.TemporaryDirectory() as t:
   sentinel=Path(t)/'important.txt';sentinel.write_text('keep')
   r=self.cli('--out',t,'--tasks','coverage');self.assertEqual(r.returncode,2);self.assertEqual(sentinel.read_text(),'keep')
 def test_invalid_concurrency_rejected(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'new';r=self.cli('--out',out,'--jobs','0','--tasks','coverage');self.assertEqual(r.returncode,2);self.assertFalse(out.exists())
 def test_all_coverage_tasks_declared(self):
  index=json.loads((R/'REPRODUCTION_INDEX.json').read_text())
  self.assertEqual(index['paper'],'v41');self.assertEqual(index['release'],'31.0.0')
  with (R/'docs/PAPER_COVERAGE.csv').open() as f:rows=list(csv.DictReader(f))
  self.assertEqual(len(rows),69)
  for row in rows:
   self.assertTrue({s.strip() for s in row['task'].split('+')}<=set(index['tasks']))
 def test_default_tasks_do_not_refit_predictors(self):
  index=json.loads((R/'REPRODUCTION_INDEX.json').read_text())['tasks']
  for name,t in index.items():
   if t.get('default',True):self.assertFalse(t.get('refits_predictors',False),name)
if __name__=='__main__':unittest.main()
