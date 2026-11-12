#!/usr/bin/env python3
"""Geometric/text PDF checks and page contact sheets for manual visual review."""
from pathlib import Path
import json,argparse
import fitz
from PIL import Image,ImageDraw
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--previews',required=True);args=ap.parse_args()
out=Path(args.previews);out.mkdir(parents=True,exist_ok=True)
reports=[]
needles=['Target-side calibration of frozen partner','Supplementary Note 21','Supplementary Table S44','Supplementary Table S45','Supplementary Table S46','Supplementary Table S47','Supplementary Table S48','Supplementary Figure S5','Comparator-inclusion sensitivity','historical pre-amendment','Functional-screen completion','Code availability','Figure 4 |','Supplementary Table S31','Supplementary Table S34','Supplementary Table S35','Supplementary Note 16','Supplementary Table S36','Supplementary reference correspondence','References','Supplementary Note 17','Supplementary Table S37','Supplementary Table S40','Supplementary Table S43']
for stem in ['WitnessCell_v36_Main_Methods','WitnessCell_v36_Supplementary_Information']:
    doc=fitz.open(R/(stem+'.pdf'));outside=[];missing=[];hits={};thumbs=[]
    for i,p in enumerate(doc):
        text=p.get_text()
        if '\ufffd' in text or '\u25a1' in text:missing.append(i+1)
        for w in p.get_text('words'):
            if w[0]<0 or w[1]<0 or w[2]>p.rect.width+1 or w[3]>p.rect.height+1:outside.append([i+1,w[4]])
        for needle in needles:
            if needle in text:hits.setdefault(needle,[]).append(i+1)
        pix=p.get_pixmap(matrix=fitz.Matrix(.55,.55),alpha=False)
        im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
        tile=Image.new('RGB',(360,490),'white');im.thumbnail((340,465));tile.paste(im,((360-im.width)//2,20))
        ImageDraw.Draw(tile).text((8,5),f'Page {i+1}',fill='black');thumbs.append(tile)
    for start in range(0,len(thumbs),12):
        contact=Image.new('RGB',(360*4,490*3),'#e6edf1')
        for j,im in enumerate(thumbs[start:start+12]):contact.paste(im,((j%4)*360,(j//4)*490))
        contact.save(out/(stem+f'_contact_{start//12+1}.png'))
    reports.append(dict(file=stem+'.pdf',pages=len(doc),words_outside_page=outside,
                        missing_glyph_markers_pages=missing,target_page_locations=hits))
    assert not outside and not missing,reports[-1]
warnings={p.name:p.read_text() for p in (R/'code').glob('WitnessCell_v36_*_render.log')}
assert not any('Missing character' in v for v in warnings.values()),warnings
audit=dict(status='PASS',documents=reports,scope='Page-bound and glyph checks; contact sheets generated for separate visual review.',
           render_warnings=warnings)
(R/'PDF_LAYOUT_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit))
