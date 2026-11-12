#!/usr/bin/env python3
"""Render readable main/Methods and SI PDFs; figures remain separate vectors."""
from pathlib import Path
import subprocess,tempfile,json,argparse,re,os
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--paper',type=Path,required=True);ap.add_argument('--only',choices=['main','si','all'],default='all');args=ap.parse_args();R=args.paper.resolve()
outputs=[]
env=os.environ.copy()
# A minimal container can contain the TeX tree but no generated format/database.
# Normal TeX installations use their existing configuration unchanged.
found=subprocess.run(['kpsewhich','xelatex.fmt'],capture_output=True,text=True)
if not found.stdout.strip() and Path('/usr/share/texlive/texmf-dist/tex/latex/base/latex.ltx').exists():
 tree='/usr/share/texlive/texmf-dist'
 runtime=tempfile.TemporaryDirectory(prefix='witnesscell_tex_')
 env.update(TEXMFCNF=tree+'/web2c',TEXMF='{'+tree+',/usr/share/texmf}',TEXINPUTS=tree+'/tex/latex/base//:'+tree+'/tex//:/usr/share/texmf/tex//:',TEXFORMATS=runtime.name+'//:',
            OPENTYPEFONTS='/usr/share/texmf/fonts/opentype//:'+tree+'/fonts/opentype//:',
            TTFONTS='/usr/share/fonts//:/usr/share/texmf/fonts/truetype//:'+tree+'/fonts/truetype//:',
            TFMFONTS='/usr/share/texmf/fonts/tfm//:'+tree+'/fonts/tfm//:',
            ENCFONTS='/usr/share/texmf/fonts/enc//:'+tree+'/fonts/enc//:',
            TEXFONTMAPS='/usr/share/texmf/fonts/map//:'+tree+'/fonts/map//:',
            TECINPUTS='/usr/share/texmf/fonts/misc//:'+tree+'/fonts/misc//:')
 fmt=subprocess.run(['xetex','-ini','-etex','-jobname=xelatex','-progname=xelatex','-interaction=nonstopmode','-output-directory='+runtime.name,'xelatex.ini'],env=env,capture_output=True,text=True)
 if fmt.returncode:raise RuntimeError(fmt.stdout[-2000:]+fmt.stderr[-2000:])
for stem in ['WitnessCell_v41_Main_Methods','WitnessCell_v41_Supplementary_Information']:
 if args.only=='main' and 'Supplementary' in stem:continue
 if args.only=='si' and 'Main_Methods' in stem:continue
 text=(R/(stem+'.md')).read_text()
 if 'Main_Methods' in stem:
  text=text.replace('## Results\n','\\Needspace{10\\baselineskip}\n\n## Results\n')
  # Keep the reference list together at readable 9-pt size.
  text=text.replace('## References\n','\\begingroup\n\\small\n\\setstretch{1.05}\n\n## References\n')+'\n\\endgroup\n'
  legends=(R/'Figure_Legends.md').read_text().replace('## Figure','\\newpage\n\n## Figure')
  legends=legends.replace('# Main figure legends\n\n\\newpage','# Main figure legends',1)
  text+='\n\n\\newpage\n\n'+legends

 else:
  text=text.replace('## WitnessCell predicts combinatorial perturbation responses through endpoint and partner learning', r'\section*{WitnessCell predicts combinatorial perturbation responses\\through endpoint and partner learning}')
  text=text.replace('## Supplementary reference correspondence\n','\\begingroup\n\\small\n\\setstretch{1.0}\n\\setlength{\\parskip}{2pt}\n\n## Supplementary reference correspondence\n')+'\n\\endgroup\n'
  text=text.replace('## Supplementary Note 13 |','\\newpage\n\n## Supplementary Note 13 |')
  text=text.replace('## Supplementary Note 16 |','\\newpage\n\n## Supplementary Note 16 |')
  text=text.replace('## Supplementary Note 21 |','\\newpage\n\n\\begingroup\n\\AtBeginEnvironment{longtable}{\\footnotesize}\n\n## Supplementary Note 21 |')
  text=text.replace('## Supplementary figure legends\n','\\endgroup\n\\newpage\n\n## Supplementary figure legends\n',1)
  text=text.replace('**B. Added prediction, primary outcome and evidence status.**','\\Needspace{24\\baselineskip}\n\n**B. Added prediction, primary outcome and evidence status.**')
  text=text.replace('### Supplementary Table S35 |','\\Needspace{20\\baselineskip}\n\n### Supplementary Table S35 |')
  for table,lines in [(47,21),(49,18),(52,24),(53,18),(54,15),(55,15),(56,15),(57,15)]:
   text=text.replace(f'### Supplementary Table S{table} |',f'\\Needspace{{{lines}\\baselineskip}}\n\n### Supplementary Table S{table} |')
  for subheading in ['**B. Expanded baseline versus original baseline.**', '**C. All seven configurations after global-reference refitting.**', '**F. Saturation-reference absolute losses.**']:
   text=text.replace(subheading,'\\Needspace{18\\baselineskip}\n\n'+subheading)
  # Keep the complete 28-task table and its final rule on the same page.
  text=re.sub(r'(\| Task \| Candidates \| PairPick[^\n]*(?:\n\|[^\n]*)+)',
      lambda m:'\\begingroup\n\\renewcommand{\\arraystretch}{0.90}\n\n'+m.group(1)+'\n\n\\endgroup',text)

 with tempfile.TemporaryDirectory() as d:
  inp=Path(d)/'input.md';inp.write_text(text)
  cmd=['pandoc',str(inp),'--from=markdown+tex_math_dollars+raw_tex','--pdf-engine=xelatex','--standalone','-V','documentclass=article','-V','fontsize=10pt','-V','geometry:margin=22mm','-V','mainfont=Liberation Serif','-V','sansfont=Liberation Sans','-V','monofont=DejaVu Sans Mono','-V','linestretch=1.12','--include-in-header='+str(HERE/'pdf_header.tex'),'--lua-filter='+str(HERE/'wrap_code.lua'),'-o',str(R/(stem+'.pdf'))]
  run=subprocess.run(cmd,text=True,capture_output=True,env=env)
  (R/(stem+'_render.log')).write_text(run.stdout+run.stderr)
  if run.returncode:raise RuntimeError(run.stderr[-4000:])
 outputs.append(stem+'.pdf');print('Rendered',stem,flush=True)
print(json.dumps({'outputs':outputs}))
