"""Reproduce the public release's final audit (Python 3 standard library).
Does not require PDF rendering dependencies. Preserves all earlier archives.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import zipfile, hashlib, json, subprocess, sys, csv, ast
from itertools import product
HERE=Path(__file__).resolve().parent
SOURCE=HERE/'Penrose_Concise_Proof_JA_v0.9_bundle.zip'
def digest(data):return hashlib.sha256(data).hexdigest()
def output(name,header,rows):
 with (HERE/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.writer(f);w.writerow(header);w.writerows(rows)
with TemporaryDirectory(prefix='penrose_release_') as tmp:
 p=Path(tmp)
 with zipfile.ZipFile(SOURCE) as z:
  assert z.testzip() is None
  manifest=json.loads(z.read('MANIFEST_SHA256.json'))
  rows=[(n,h,digest(z.read(n))) for n,h in sorted(manifest.items())]
  assert all(h==g for n,h,g in rows)
  z.extractall(p)
 output('v09_preservation.csv',['file','expected_sha256','actual_sha256','status'],[(n,h,g,'PASS') for n,h,g in rows])
 subprocess.run([sys.executable,str(p/'review_v09.py')],cwd=p,check=True)
 with zipfile.ZipFile(SOURCE) as z:
  expected=list(csv.DictReader(z.read('replay_SHA256.csv').decode('utf-8-sig').splitlines()))
  reread=[(r['name'],digest(z.read(r['name'])),digest((p/r['name']).read_bytes())) for r in expected]
  assert all(h==g for _,h,g in reread)
 output('v09_reexecution.csv',['file','original_sha256','new_sha256','status'],[(n,h,g,'PASS') for n,h,g in reread])
 sys.path.insert(0,str(p/'_original_v14'))
 import independent as M
 states,arcs,pruned=M.independent_graph();assert len(states)==285
 # Check rank and domain without symmetry quotient; no finite-level assumption at this stage.
 # A child of depth d>0 uses a smaller index. For maximum m>=3 all formal children
 # have positive generation; K1, if ever produced, is empty and contributes nothing.
 def check_pair(X,Y):
  a,i,*_=X;b,j,*_=Y
  delta=min(i,j);assert delta in (1,2)
  st=M.pair(X,Y); assert min(st[1],st[3])==0 and max(st[1],st[3])<=1
  assert (st[1],st[3])==(i-delta,j-delta)
  assert 3-delta>=1 and 3-i>=1 and 3-j>=1
  return delta,st
 rows=[];total=0
 for st in sorted(states):
  a,i,b,j,e,r,d=st
  assert min(i,j)==0 and max(i,j)<=1
  X=M.desc(a,i);Y=M.desc(b,j,e,r,d)
  xs=M.expand(X) if i==0 else [X];ys=M.expand(Y) if j==0 else [Y]
  counts={1:0,2:0};outside=0
  for U,V in product(xs,ys):
   delta,ch=check_pair(U,V);counts[delta]+=1
   if M.aoutside(ch[0],ch[1],ch[2],ch[3],ch[-1]):outside+=1
   else:assert ch in states
  total+=len(xs)*len(ys)
  rows.append((repr(st),counts[1],counts[2],outside,'PASS'))
 assert total==arcs==20855 and sum(r[3] for r in rows)==pruned
 output('rank_domain_285_states.csv',['state','drop1','drop2','disk_prune','status'],rows)
 # Parent grammar roots: guarantee origin and child rank, including cross-generation roots.
 from itertools import combinations
 rootrows=[]
 for t in 'ABKL':
  for (i,X),(j,Y) in combinations(enumerate(M.expand(M.desc(t))),2):
   delta,ch=check_pair(X,Y)
   out=M.aoutside(ch[0],ch[1],ch[2],ch[3],ch[-1]);assert out or ch in states
   rootrows.append((t,i,j,delta,repr(ch),'DISK_PRUNE' if out else 'STATE','PASS'))
 assert len(rootrows)==130
 output('rank_domain_130_roots.csv',['parent','child1','child2','drop','state','result','status'],rootrows)
 base=M.BASE
 assert base['K',1]==set() and base['K',2]=={((-1,0,0,0),1)} and base['L',1]=={((-1,0,0,0),-1),((0,-1,0,0),-1),((0,0,0,0),1),((1,1,1,1),-1)}
 audit=json.loads((p/'audit_results.json').read_text())
 audit.update(v09_manifest=len(manifest),v09_replayed_files=len(reread),rank_states=len(states),rank_transitions=total,rank_drop1=sum(r[1] for r in rows),rank_drop2=sum(r[2] for r in rows),rank_roots=len(rootrows),source_sha256=digest(SOURCE.read_bytes()))
 (HERE/'final_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
 print(audit)
