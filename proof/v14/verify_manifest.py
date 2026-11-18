from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parent
m=json.loads((r/'MANIFEST_SHA256.json').read_text())
for f,h in m.items():assert hashlib.sha256((r/f).read_bytes()).hexdigest()==h,f
print('manifest',len(m),'PASS')
