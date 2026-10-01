from pathlib import Path
import zipfile,json,re,shutil
h=Path(__file__).resolve().parent;d=h/'manuscript'
(h/'README_zh.md').write_text('''# v6阅读入口\n\n先读framework_zh.md，再读manuscript/manuscript_v6.pdf（26页，六主图、四补图）；DOCX为同源可编辑稿。host/results_zh.md记录本轮实际完成与边界。\n\n本包包括S1–S3证据表、当前图件和候选源数据子集，供科学讨论及改稿使用。不是完整历史SI、完整分析复现环境或最终投稿包。原45–47群落读图不并入，原v5及学生稿未覆盖。\n\n最终页面来自render4；本目录早期render、render2、render3仅保留内部版本历史，未收入阅读包。新增候选图沿用原数据与几何，只纠正标签。\n\n尚待整合：完整历史SI、作者缺失单位2、公共代码/数据归档、目标期刊格式。生态机制与新蛋白功能均未被本轮确证。\n''')
paths=[h/'framework_zh.md',h/'README_zh.md',h/'host/results_zh.md',h/'host/delivery_qa.json',h/'independent_audit/dpfq_evidence.md',h/'independent_audit/review_v6.md']
paths += [d/f'manuscript_v6.{ext}' for ext in ['md','docx','pdf']]
paths += list(d.glob('Supplementary*.tsv'))+list((d/'Supplementary_Table_S3_source_data').glob('*'))
s=(d/'manuscript_v6.md').read_text();paths += [d/r for r in re.findall(r'!\[.*?\]\((assets/[^)]+)\)',s)]
paths += list((h/'figures').glob('figure6*'))
with zipfile.ZipFile(h/'reading_bundle_v6.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in paths:assert p.is_file(),p;z.write(p,p.relative_to(h))
with zipfile.ZipFile(h/'reading_bundle_v6.zip') as z:assert z.testzip() is None;print('reading bundle',len(z.namelist()),'files')
