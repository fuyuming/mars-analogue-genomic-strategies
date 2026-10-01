from pathlib import Path
import json,re,zipfile,hashlib,shutil,xml.etree.ElementTree as ET
from collections import Counter
import fitz
h=Path(__file__).resolve().parent;d=h/'manuscript';prior=h.parent/'51_evidence_figures_20261001/manuscript'
shutil.copy2(d/'render4/manuscript_v6.pdf',d/'manuscript_v6.pdf')
s=(d/'manuscript_v6.md').read_text();ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with zipfile.ZipFile(d/'manuscript_v6.docx') as z:
 root=ET.fromstring(z.read('word/document.xml'))
 actual=[''.join(p.itertext()) for p in []]
 actual=[''.join(t.text or '' for t in p.findall('.//w:t',ns)) for p in root.findall('./w:body/w:p',ns)]
 actual=[p for p in actual if p]
 imgs=[n for n in z.namelist() if n.startswith('word/media/')]
 links=[e.attrib.get('{'+ns['w']+'}anchor') for e in root.findall('.//w:hyperlink',ns) if e.attrib.get('{'+ns['w']+'}anchor')]
 bookmarks={e.attrib.get('{'+ns['w']+'}name') for e in root.findall('.//w:bookmarkStart',ns)}
expected=[]
for b in s.split('\n\n'):
 b=b.strip()
 if not b or b=='## Figures and legends' or b.startswith('!['): continue
 expected.append(re.sub('^#{1,3} ','',b).replace('\n',' '))
assert Counter(expected)==Counter(actual), {'missing':list((Counter(expected)-Counter(actual)).elements()),'extra':list((Counter(actual)-Counter(expected)).elements())}
assert len(imgs)==10 and set(links)<=bookmarks
assert len(re.findall(r'^\d+\. ',s,re.M))==43
unchanged=[]
for rel in re.findall(r'!\[.*?\]\((assets/[^)]+)\)',s):
 if 'figure6_' not in rel:
  assert (d/rel).read_bytes()==(prior/rel).read_bytes(),rel;unchanged.append(rel)
for p in d.glob('Supplementary_Table_S[12]*tsv'):assert p.read_bytes()==(prior/p.name).read_bytes(),p
pdf=fitz.open(d/'manuscript_v6.pdf');assert len(pdf)==26
text='\n'.join(p.get_text() for p in pdf)
for term in ['0.815','0.114','648','97–230','250–359','59–69','22,369','Supplementary Table S3']:
 assert term in text,term
with (d/'Supplementary_Table_S3_candidate_evidence.tsv').open() as f:
 import csv
 rows=list(csv.DictReader(f,delimiter='\t'))
for row in rows:assert (d/'Supplementary_Table_S3_source_data'/row['retained_source']).exists(),row
bounds=[]
for i,page in enumerate(pdf):
 for b in page.get_text('blocks'):
  if b[0]<-1 or b[1]<-1 or b[2]>page.rect.width+1 or b[3]>page.rect.height+1: bounds.append(i+1)
assert not bounds,bounds
r={'MD_DOCX_paragraph_multiset_exact':True,'DOCX_images':len(imgs),'PDF_pages':len(pdf),'reference_entries':43,'native_citation_targets_resolve':True,'all_26_pages_host_visual_contact_review':True,'text_outside_page_bounds':bounds,'nine_inherited_figures_byte_unchanged':unchanged,'S1_S2_byte_unchanged':True,'S3_index_rows':len(rows),'final_render':'render4','scope':'Existing-evidence integration, no new biological analysis or full historical pipeline replication; reading draft, not complete submission package','file_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [d/'manuscript_v6.md',d/'manuscript_v6.docx',d/'manuscript_v6.pdf']}}
(h/'host/delivery_qa.json').write_text(json.dumps(r,ensure_ascii=False,indent=2));print(json.dumps(r,ensure_ascii=False,indent=2))
