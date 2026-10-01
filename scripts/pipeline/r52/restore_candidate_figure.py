from pathlib import Path
import shutil,subprocess,xml.etree.ElementTree as E,json,hashlib
H=Path(__file__).resolve().parent;A=H.parents[1].parent/'code/result_raw/ISME_Figure3_DPFQ008_lineage_context_20260729_v11';F=H/'figures';D=H/'manuscript'
# Edit original vector text only. Geometry, data, tree and plotted statistics unchanged.
p=A/'Figure3_DPFQ008_lineage_context.svg';s=p.read_text()
changes={'cytochrome c-like ×2':'cytochrome c-like','Non-random genomic context':'Context in selected carriers','A cross-environment but lineage-structured family':'A lineage-structured candidate across sources','non-Qaidam pairs:':'pairs excluding Q–Q:'}
for a,b in changes.items():assert a in s,a;s=s.replace(a,b)
# Remove forced text stretching for edited labels only.
root=E.fromstring(s)
for el in root.iter():
 if el.tag.endswith('text') and any(b in (el.text or '') for b in changes.values()):
  el.attrib.pop('textLength',None);el.attrib.pop('lengthAdjust',None)
svg=F/'figure6_DPFQ008.svg';E.register_namespace('','http://www.w3.org/2000/svg');E.ElementTree(root).write(svg,encoding='utf-8',xml_declaration=True)
subprocess.run(['/usr/local/bin/rsvg-convert','-f','pdf','-o',str(F/'figure6_DPFQ008.pdf'),str(svg)],check=True)
subprocess.run(['/usr/local/bin/rsvg-convert','-w','2600','-o',str(F/'figure6_DPFQ008.png'),str(svg)],check=True)
shutil.copy2(F/'figure6_DPFQ008.png',D/'assets/figure6_DPFQ008.png')
S=D/'Supplementary_Table_S3_source_data';S.mkdir(exist_ok=True)
for f in (A/'Source_Data').glob('*'):shutil.copy2(f,S/f.name)
(H/'host/candidate_figure_provenance.json').write_text(json.dumps({'original_svg':str(p.resolve()),'original_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'changes':changes,'no_data_or_geometry_change':True},indent=2,ensure_ascii=False))
