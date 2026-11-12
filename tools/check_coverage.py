from pathlib import Path
import re,csv,json,argparse,hashlib
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.resolve().exists() and not a.out.resolve().is_relative_to(R);a.out.mkdir(parents=True)
rows=list(csv.DictReader((R/'docs/PAPER_COVERAGE.csv').open()));tables={int(x['item'][1:]):x for x in rows if re.fullmatch(r'S\d+',x['item'])};si=(R/'paper/manuscript_sources/WitnessCell_v41_Supplementary_Information.md').read_text();titles={int(n):t for n,t in re.findall(r'^### Supplementary Table S(\d+) \| (.+)$',si,re.M)};assert set(tables)==set(titles)==set(range(1,58))
for n in tables:assert tables[n]['title']==titles[n]
for row in rows:
 for item in row['source'].split(' + '):assert (R/item).exists(),item
 assert all(x.strip() in json.loads((R/'REPRODUCTION_INDEX.json').read_text())['tasks'] for x in row['task'].split('+'))
fig_text=(R/'paper/manuscript_sources/Figure_Legends.md').read_text();sup_text=(R/'paper/manuscript_sources/Supplementary_Figure_Legends.md').read_text()
fig_titles={f'Fig{n}':v for n,v in re.findall(r'^## Figure (\d+) \| (.+)$',fig_text,re.M)}
fig_titles.update({f'SuppFig{n}':v for n,v in re.findall(r'^## Supplementary Figure S(\d+) \| (.+)$',sup_text,re.M)})
assert len(fig_titles)==12
for row in rows:
 if row['item'] in fig_titles:assert row['title']==fig_titles[row['item']]
assert len(rows)==69
report={'status':'PASS','supplementary_tables':57,'main_figures':5,'supplementary_figures':7,'mapped_items':len(rows),'scope':'Complete current item/title/source/task mapping; numerical verification is performed by the mapped task, not inferred from presence.'}
(a.out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
