#!/usr/bin/env python3
"""Count manuscript prose without dropping the first introduction paragraph."""
from pathlib import Path
import json,re
R=Path(__file__).resolve().parents[1]

def counts(text):
    # Abstract is one paragraph; the unheaded introduction starts immediately after it.
    body=text.split('## Abstract\n',1)[1].lstrip('\n')
    abstract,main=body.split('\n\n',1)
    main=main.split('## Methods\n',1)[0]
    headings_removed=re.sub(r'^#{1,6}[^\n]*\n?', '', main, flags=re.M)
    sections={}
    for name,chunk in [('Introduction',main.split('## Results\n')[0]),
                       ('Results',main.split('## Results\n')[1].split('## Discussion\n')[0]),
                       ('Discussion',main.split('## Discussion\n')[1])]:
        sections[name]=len(re.sub(r'^#{1,6}[^\n]*\n?', '',chunk,flags=re.M).split())
    references=text.split('## References\n',1)[1]
    return dict(abstract_words=len(abstract.split()),
                main_prose_whitespace_words=len(headings_removed.split()),
                sections=sections,
                counting_rule='Whitespace tokens; title, abstract, section headings, Methods, references and figure legends excluded from main prose; displayed equations retained. The first unheaded introduction paragraph is included.',
                references=len(re.findall(r'^\d+\. ',references,re.M)))

if __name__=='__main__':
    old=counts((R/'code/source_inputs/WitnessCell_v26_Main_Methods.md').read_text())
    assert old['main_prose_whitespace_words']==4155 and old['abstract_words']==128
    text=(R/'WitnessCell_v36_Main_Methods.md').read_text()
    report=counts(text)
    # A sentinel in the first unheaded paragraph must affect prose, not the abstract.
    changed=text.replace('Combinatorial perturbation screens','COUNT_SENTINEL Combinatorial perturbation screens',1)
    delta=counts(changed)
    assert delta['main_prose_whitespace_words']==report['main_prose_whitespace_words']+1
    assert delta['abstract_words']==report['abstract_words']
    assert sum(report['sections'].values())==report['main_prose_whitespace_words']
    assert report['abstract_words']<=150 and report['main_prose_whitespace_words']<=5000
    (R/'WORD_COUNTS.json').write_text(json.dumps(report,indent=2)+'\n')
    (R/'revision/WORD_COUNT_REGRESSION.json').write_text(json.dumps(dict(status='PASS',
        v26_corrected_prose=4155,v26_previously_reported_prose=4064,
        first_introduction_paragraph_sentinel_counted=True,section_sum_matches=True),indent=2)+'\n')
    print(json.dumps(report))
