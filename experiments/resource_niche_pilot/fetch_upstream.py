"""Fetch pinned upstream code/data without importing or running the program."""
from pathlib import Path
import hashlib,json,urllib.request
B=Path(__file__).resolve().parent
for row in json.loads((B/'upstream_manifest.json').read_text()):
    p=B/row['path'];p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():
        data=p.read_bytes()
    else:
        with urllib.request.urlopen(row['url'],timeout=60) as response:
            data=response.read()
    assert len(data)==row['bytes'],row['path']
    assert hashlib.sha256(data).hexdigest()==row['sha256'],row['path']
    if not p.exists():p.write_bytes(data)
    print('Verified',row['path'])
