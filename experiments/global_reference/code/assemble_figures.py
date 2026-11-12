"""Assemble current vector panels without rasterizing them."""
from pathlib import Path
import argparse, re
from pypdf import PdfReader, PdfWriter

ap = argparse.ArgumentParser(); ap.add_argument('--paper', type=Path, required=True); a = ap.parse_args()
for folder, label, expected in [('figures', 'Main', 5), ('supplementary_figures', 'Supplementary', 7)]:
    files = sorted((a.paper / folder).glob('*.pdf'), key=lambda p: int(re.search(r'(?:Supp)?Fig(\d+)', p.name)[1]))
    assert len(files) == expected
    writer = PdfWriter()
    for source in files:
        reader = PdfReader(source); assert len(reader.pages) == 1
        writer.add_page(reader.pages[0])
    writer.add_metadata({'/Title': f'WitnessCell v39 {label} Figures'})
    with (a.paper / f'WitnessCell_v39_{label}_Figures.pdf').open('wb') as f:
        writer.write(f)
    print(label, expected)
