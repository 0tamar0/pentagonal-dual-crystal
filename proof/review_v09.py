"""Audit from the preserved v0.4 archive. Independent polygon test, shared exact generation.
Python 3 standard library; no old proof outputs taken as an assertion of validity.
"""
import ast,csv,hashlib,io,json,sys,zipfile,math
from pathlib import Path
from itertools import combinations,product
from functools import lru_cache
from collections import Counter
O=Path(__file__).resolve().parent
src=Path(sys.argv[1]) if len(sys.argv)>1 else O/'Penrose_Concise_Proof_JA_v0.4_bundle.zip'
def unpack(data):
 z=zipfile.ZipFile(io.BytesIO(data));assert z.testzip() is None
 m=json.loads(z.read('MANIFEST_SHA256.json'))
 assert all(hashlib.sha256(z.read(n)).hexdigest()==h for n,h in m.items())
 return z,len(m)
v4,k4=unpack(src.read_bytes());v3,k3=unpack(v4.read('Penrose_Concise_Proof_JA_v0.3_bundle.zip'));v14,k14=unpack(v3.read('Penrose_AaBb_Induction_1000Go_v14_20261119.zip'))
P=O/'_original_v14';P.mkdir(exist_ok=True);v14.extractall(P);sys.path.insert(0,str(P))
from independent import *
def read(z,n):return list(csv.DictReader(io.StringIO(z.read(n).decode('utf-8-sig'))))
def write(n,head,rows):
 with (O/n).open('w',encoding='utf-8-sig',newline='') as f:w=csv.writer(f);w.writerow(head);w.writerows(rows)
# Check exact specified grammar, containment and symmetry, without executing old regression scripts.
assert len([x for t in 'ABKL' for x in GR[t]])==32
for (t,n),S in BASE.items():
 assert {(conj(c),e) for c,e in S}==S
 if t in 'AB':assert all(mapset(S,Z,1,j)==S for j in range(5))
for t,g in GR.items():
 for u,d,e,j,c in g:
  assert any(sym_eq((u,d,e,-j%5,conj(c)),Y) for Y in g)
  if t in 'AB':
   for r in range(5):assert any(sym_eq((u,d,e,(j+r)%5,rot(c,1,r)),Y) for Y in g)
subs=[r for t in 'KL' for r in subset(t)];cov=[r for t in 'KL' for r in coverage(t)]
assert subs==[tuple([a['gap'],int(a['retained_child']),a['child'],int(a['source_child']),a['source']]) for a in read(v14,'subset_10_identities.csv')]
assert cov==[tuple([a['gap'],int(a['source_child']),a['source'],a['case'],a['witness_path'],a['target']]) for a in read(v14,'coverage_22_identities.csv')]
for n in range(1,5):
 for t in 'ABKL':
  G=actual(t,n)
  if t in 'AB':R=RAW[t,n]
  elif t=='L':R=RAW['A',n]-trans(RAW['B',n],mul(neg(pp(n)),inv),-1)
  else:R=set() if n==1 else RAW['B',n-1]-trans(RAW['A',n],neg(pp(n-1)),1)
  assert G==R
# Rebuild literal unquotiented graph and every canonical transition, not just totals.
seen,arcs,pruned=independent_graph();old={ast.literal_eval(a['state']) for a in read(v14,'independent_unquotiented_states.csv')};assert seen==old
roots,states,ar,pr=abuild();sid={st:i for i,st in enumerate(states)}
assert {acanon(*st) for st in seen}==set(states)
oldtr=Counter((int(a['parent']),int(a['level_drop']),ast.literal_eval(a['child']),a['destination_or_prune']) for a in read(v14,'atom_transition_certificate.csv'))
newtr=Counter((sid[st],d,ch,'DISK_PRUNE' if aoutside(*ch[:4],ch[-1]) else str(sid[ch])) for st in states for ch,d in achildren(st));assert oldtr==newtr
print('manifests',k4,k3,k14,'graph',len(seen),arcs,pruned,'inclusion',len(subs),len(cov),flush=True)
# New polygon kernel: cyclic convolution of five integer coefficients. Directed edges,
# cross-product signs, and exact tangent segment endpoints. No classify()/interval() used.
def c5(a):return (*a,0)
def sub5(a,b):return tuple(x-y for x,y in zip(a,b))
def norm5(a):return tuple(x-a[4] for x in a)
@lru_cache(None)
def prod5(a,b):return norm5(tuple(sum(a[i]*b[(k-i)%5] for i in range(5)) for k in range(5)))
def conj5(a):return (a[0],a[4],a[3],a[2],a[1])
def sgn(a,b):
 x=2*a+b;y=b
 if x==0:return (y>0)-(y<0)
 if y==0 or (x>0)==(y>0):return (x>0)-(x<0)
 return ((x>0)-(x<0)) if x*x>5*y*y else ((y>0)-(y<0))
def cross5(a,b):
 c=prod5(conj5(a),b);v=c[2]-c[3];return sgn(c[1]-c[4]-v,v)
def dot5(a,b):
 c=prod5(conj5(a),b);v=c[1]+c[4];return sgn(2*c[0]-v,v-c[2]-c[3])
@lru_cache(None)
def vertices5(T):return tuple(c5(x) for x in verts(T))
@lru_cache(None)
def rel(d,e,f):
 if d==Z and e==f:return 'identical'
 # Quick exact separation by radius. This is independent of old polygon classifier.
 if dot5(c5(d),c5(d)) and dot5(c5(d),c5(d))>0:
  # |d|^2 > (2/phi)^2 using five-coefficient multiplication
  r=prod5(c5(d),conj5(c5(d)));lim=prod5(c5(times(inv,2)),c5(times(inv,2)))
  if dot5(sub5(r,lim),c5(I))>0:return 'separate'
 X=vertices5((Z,e));Y=vertices5((d,f));tangent=None
 for P,Q in ((X,Y),(Y,X)):
  for v,w in zip(P,P[1:]+P[:1]):
   a=sub5(w,v);signs=[cross5(a,sub5(x,v)) for x in Q]
   if max(signs)<0:return 'separate'
   if max(signs)==0:tangent=v,w,tuple(x for x,k in zip(Q,signs) if k==0)
 if tangent is None:return 'BAD_area'
 v,w,T=tangent;a=sub5(w,v)
 def cmp(x,y):return dot5(sub5(x,y),a)
 qlo,qhi=T[0],T[-1]
 if cmp(qlo,qhi)>0:qlo,qhi=qhi,qlo
 lo=qlo if cmp(qlo,v)>0 else v;hi=qhi if cmp(qhi,w)<0 else w;k=cmp(hi,lo)
 if k<0:return 'separate'
 if k==0:return 'point' if lo in (v,w) and lo in T else 'BAD_partial_point'
 return 'edge' if lo==v==qlo and hi==w==qhi and len(T)==2 else 'BAD_partial_edge'
def relation(x,y):return rel(sub(y[0],x[0]),x[1],y[1])
# Execute all unsymmetrized terminal cases. Compare each exact classification with original.
cat=set();rr=[]
for st in sorted(seen):
 a,i,b,j,e,r,d=st
 for n in (1,2):
  if n<=max(i,j):continue
  X=atoms(a,n-i);Y=[(t,mapset(T,mul(d,pp(n-1)),e,r)) for t,T in atoms(b,n-j)];C=Counter();contacts=0
  for (u,S),(v,T) in product(X,Y):
   c=Counter(relation(x,y) for x,y in product(S,T));assert not any(k=='identical' or k.startswith('BAD') for k in c),(st,n,c)
   C.update(c)
   if c['edge'] or c['point']:cat.add(signature(u,S,v,T));contacts+=1
  rr.append((str(st),n,sum(C.values()),C['edge'],C['point'],contacts))
oldcat={ast.literal_eval(a['full_two_atom_coordinate_signature']):a['type'] for a in read(v14,'10_full_exact_canonical_motifs.csv')};assert cat==set(oldcat) and len(rr)==517
write('517_polygon_rechecks.csv',['state','m','tile_pairs','edge_pairs','point_pairs','atom_contacts'],rr)
print('polygon audit',len(rr),sum(r[2] for r in rr),len(cat),flush=True)
# Verify original small assemblies, including atomic interiors, and motif realization in A3 union B3.
occ=[]
for t,n in BASE:
 T=BASE[t,n];assert all(relation(x,y) in ('separate','point','edge') for x,y in combinations(T,2))
for n in (1,2,3):
 for t in 'ABKL':
  A=atoms(t,n);C=Counter();checks=0
  for (u,S),(v,T) in combinations(A,2):
   c=Counter(relation(x,y) for x,y in product(S,T));checks+=sum(c.values());assert not any(k=='identical' or k.startswith('BAD') for k in c)
   if c['edge'] or c['point']:C[oldcat[signature(u,S,v,T)]]+=1
  occ.append((t,n,len(A),checks,str(sorted(C.items())),sum(C.values())))
write('12_assembly_rechecks.csv',['type','n','atoms','tile_pairs','motifs','contact_pairs'],occ)
# Validate 662 diagram child pairs and full sets from generating coordinates.
@lru_cache(None)
def incidence(T):
 V=set();E=set()
 for x in T:
  vv=vertices5(x);V.update(vv);E.update(tuple(sorted((a,b))) for a,b in zip(vv,vv[1:]+vv[:1]))
 return E,V
parts={}
for a in read(v4,'figure_parts.csv'):
 t={'a':'K','b':'L'}.get(a['type'],a['type']);n=int(a['generation']);e=int(a['sign']);r=int(a['rotation_72deg']);c=ast.literal_eval(a['translation_normalized']);T=mapset(actual(t,n),mul(c,pp(2)),e,r);assert len(T)==int(a['tiles']);parts[a['rule'],a['id']]=frozenset(T)
fr=[]
for a in read(v4,'figure_child_pairs.csv'):
 X,V=incidence(parts[a['rule'],a['left_id']]);Y,W=incidence(parts[a['rule'],a['right_id']]);E=X&Y;P=V&W;got=(len(E),len(P-{x for edge in E for x in edge}),len(P));assert got==tuple(int(a[x]) for x in ('shared_edges','isolated_points','shared_vertices'));fr.append((a['rule'],a['left_id'],a['right_id'],*got,'PASS'))
write('662_figure_pairs.csv',['rule','left','right','edges','isolated_vertices','vertices','status'],fr)
# Check 105,990 complete finite base pairs (including within atoms) and three genuine bad intersections.
bas=[]
for n in (1,2,3):
 for t in 'AB':
  S=RAW[t,n];c=Counter(relation(x,y) for x,y in combinations(S,2));assert not any(k=='identical' or k.startswith('BAD') for k in c);bas.append((t,n,len(S),sum(c.values()),c['separate'],c['edge'],c['point']))
assert sum(x[3] for x in bas)==105990
write('105990_base_polygon_pairs.csv',['type','n','tiles','pairs','separate','edge','point'],bas)
negatives=[]
for a in read(v14,'explicit_geometric_negative_controls.csv'):
 ret=relation(ast.literal_eval(a['test_polygon1']),ast.literal_eval(a['test_polygon2']));assert ret==a['classification'];negatives.append((a['case'],ret,'REJECTED'))
write('negative_controls.csv',['case','result','status'],negatives)
# Re-evaluate activity, all 14 standard states and the asymmetric B-b grouping.
active={st for st in states if any(abase(st,n,set()) for n in (1,2))}
while True:
 nxt=active|{st for st in states if any(ch in active for ch,d in achildren(st))}
 if nxt==active:break
 active=nxt
assert len(active)==13 and sum(st in active for t,i,j,k,st in roots)==44
rules={a['rule']:ast.literal_eval(a['state']) for a in read(v14,'14_standard_and_auxiliary_rules.csv')};names={v:k for k,v in rules.items()};rrules=[]
for a in read(v14,'14_standard_and_auxiliary_rules.csv'):
 st=rules[a['rule']];ch=list(achildren(st));assert all(u in sid or aoutside(*u[:4],u[-1]) for u,d in ch)
 ct=Counter(names[u] for u,d in ch if u in active);assert ct==Counter(dict(ast.literal_eval(a['nonzero_child_rules'])))
 for n in (1,2):
  x,i,y,j,e,r,d=st
  if n<=max(i,j):continue
  S=actual(x,n-i);T=mapset(actual(y,n-j),mul(d,pp(n-1)),e,r);cs=Counter(relation(v,w) for v,w in product(S,T));assert not any(k.startswith('BAD') or k=='identical' for k in cs)
 rrules.append((a['rule'],len(ch),sum(ct.values()),str(sorted(ct.items())),'PASS'))
write('14_rule_rechecks.csv',['rule','all_candidates','active','child_rules','status'],rrules)
for a in read(v14,'Bb_2_plus_1_exact_grouping.csv'):
 st=rules[a['rule']];x,i,y,j,e,r,d=st;U=desc(x,i);W=expand(desc(y,j,e,r,d))[int(a['gap_child'])];q=acanon(*pair(U,W));assert q==ast.literal_eval(a['pair_state'])
 ct=Counter(names[c] for c,k in achildren(q) if c in active);assert ct==Counter(dict(ast.literal_eval(a['expanded_active_children'])))
# Recompute each SVG's polygon vertices in display coordinates. Float is only for figure layout.
import cmath,xml.etree.ElementTree as ET
zz0=cmath.exp(2j*math.pi/5)
def xy(q):return sum(v*zz0**i for i,v in enumerate(q))
svgrows=[]
for rule,st in rules.items():
 l,i,r,j,e,k,d=st;X=actual(l,3-i);Y=mapset(actual(r,3-j),mul(d,pp(2)),e,k)
 P=[S for (name,ident),S in parts.items() if name==rule and ident.startswith('X')];Q=[S for (name,ident),S in parts.items() if name==rule and ident.startswith('Y')]
 assert set().union(*P)==X and set().union(*Q)==Y and sum(map(len,P))==len(X) and sum(map(len,Q))==len(Y)
 vs=[xy(v) for t in X|Y for v in verts(t)];xmin,xmax=min(v.real for v in vs),max(v.real for v in vs);ymin,ymax=min(v.imag for v in vs),max(v.imag for v in vs);sc=min(392/(xmax-xmin),248/(ymax-ymin));mx=(xmin+xmax)/2;my=(ymin+ymax)/2
 def pt(q):
  c=xy(q);return (210+(c.real-mx)*sc,138-(c.imag-my)*sc)
 def poly(t):return tuple((round(x,3),round(y,3)) for v in verts(t) for x,y in [pt(v)])
 expected=Counter(poly(t) for t in X|Y)
 for mode in ('before','after'):
  root=ET.fromstring(v4.read('figures/'+rule+'_'+mode+'.svg'));got=Counter(tuple(tuple(map(float,p.split(','))) for p in a.attrib['points'].split()) for a in root.findall('{http://www.w3.org/2000/svg}polygon'));assert expected==got,(rule,mode)
  svgrows.append((rule,mode,sum(got.values()),'PASS'))
write('28_svg_coordinate_rechecks.csv',['rule','panel','polygons','status'],svgrows)

summary=dict(manifests=[k4,k3,k14],inclusion=len(subs)+len(cov),unquotiented_states=len(seen),arcs=arcs,pruned=pruned,canonical_states=len(states),canonical_arcs=ar,base_cases=len(rr),base_tile_pairs=sum(r[2] for r in rr),motifs=len(cat),figure_pairs=len(fr),figure_parts=len(parts),svg_panels=len(svgrows),svg_polygons=sum(r[2] for r in svgrows),negative_controls=len(negatives),standard_rules=len(rrules),base_full_pairs=sum(r[3] for r in bas))
(O/'audit_results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(summary)
