"""Reproduce the target-only extraction from the existing original KO calls.
No gene prediction or annotation is run; output is refused if it already exists.
"""
from pathlib import Path
import argparse,csv,gzip,json
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
h=Path(__file__).resolve().parent
rules=json.loads((h.parent/'49_system_configuration_20261001/host/system_rules_provisional.json').read_text())['rules']
keep={k for ks in rules.values() for k in ks};seen=selected=0
with gzip.open(a.input,'rt') as inp,a.output.open('x') as out:
 reader=csv.DictReader(inp,delimiter='\t');writer=csv.DictWriter(out,fieldnames=reader.fieldnames,delimiter='\t');writer.writeheader()
 for row in reader:
  seen+=1
  if row['KO'] in keep:writer.writerow(row);selected+=1
print(json.dumps({'input_rows':seen,'target_rows':selected}))
