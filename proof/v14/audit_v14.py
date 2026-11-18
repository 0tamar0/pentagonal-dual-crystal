"""Generate exact graph, coverage, base, catalog and 1000-regression certificates."""
import hashlib,zipfile,io,json
from independent import *
from collections import Counter

def write(name,header,rows):
 with (ROOT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.writer(f);w.writerow(header);w.writerows(rows)
def read(name):
 with (ROOT/name).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def oldhash(sig):
 a,b,X,Y=sig; labels={'A':'A1','B':'B1','K':'a2','L':'b1'}
 a,b=labels[a],labels[b]
 # canonicalize again with historical orientation encoding (0=+,1=-) and label order.
 if a>b:a,b,X,Y=b,a,Y,X
 options=[]
 for mirror in (0,1):
  for e in (1,-1):
   for j in range(5):
    def F(S):return sorted((rot(conj(p) if mirror else p,e,j),int(f*e==-1)) for p,f in S)
    M,N=F(X),F(Y);o=min(p for p,f in M+N)
    U=tuple((sub(p,o),f) for p,f in M);V=tuple((sub(p,o),f) for p,f in N)
    if a==b and U>V:U,V=V,U
    options.append((a,b,U,V))
 return hashlib.sha256(repr(min(options)).encode()).hexdigest()

# Verify original input archive in place, not only its display name.
prior=ROOT/'Penrose_AaBb_Induction_100Go_v13_20261119.zip'
with zipfile.ZipFile(prior) as zf:
 assert zf.testzip() is None
 manifest=json.loads(zf.read('MANIFEST_SHA256.json'))
 # v13 uses filename -> sha256; allow only its actual format.
 priorchecks=[]
 for f,h in manifest.items():
  got=hashlib.sha256(zf.read(f)).hexdigest();assert got==h,(f,got,h)
  priorchecks.append((f,h,got,'PASS'))
 assert zf.read('10_base_motifs.csv')==(ROOT/'prior_10_motifs.csv').read_bytes()
write('prior_v13_manifest_check.csv',['name','expected','actual','status'],priorchecks)

# Exact regular-pentagon bases and symmetry hypotheses.
base_rows=[]
for (t,n),S in sorted(BASE.items()):
 assert mapset(S,Z,1,0)==S
 assert {(conj(p),e) for p,e in S}==S
 if t in 'AB':assert all(mapset(S,Z,1,j)==S for j in range(5))
 rel=Counter(classify(x,y) for x,y in combinations(S,2));assert not any(k for u,k in rel.items() if u.startswith('BAD') or u=='identical')
 base_rows.append((t,n,len(S),sum(rel.values()),rel['separate'],rel['point'],rel['edge'],'PASS'))
# Reflection closure of child descriptors, and rotational closure for A/B.
for t,G in GR.items():
 for u,d,e,j,c in G:
  X=(u,d,e,-j%5,conj(c));assert any(sym_eq(X,Y) for Y in G)
  if t in 'AB':
   for r in range(5):assert any(sym_eq((u,d,e,(j+r)%5,rot(c,1,r)),Y) for Y in G)
write('atomic_base_SAT.csv',['kind','n','tiles','pairs','separate','point','edge','status'],base_rows)
# The raw 15-state graph with independent base SAT (allows intentionally identical copies in unions).
rroots,rstates,rarcs,rpruned=build();rsid={st:i for i,st in enumerate(rstates)};rrows=[];rbase=[];rt=[]
BAS=gen(2)
for st in rstates:
 i=rsid[st];rrows.append((i,*st))
 for ch,drop in children(st):
  a,ii,b,j,e,d=ch;outside_=outside(a,ii,b,j,d)
  assert outside_ or ch in rsid;rt.append((i,drop,str(ch),'DISK_PRUNE' if outside_ else rsid[ch]))
 for n in (1,2):
  rel=base(st,n,BAS)
  if rel is None:continue
  assert not any(k for u,k in rel.items() if u.startswith('BAD'))
  rbase.append((i,n,sum(rel.values()),rel['identical'],rel['separate'],rel['point'],rel['edge'],'PASS'))
write('raw_15_states.csv',['id','a','depth_a','b','depth_b','sign','d'],rrows)
write('raw_transition_certificate.csv',['parent','level_drop','child','destination_or_prune'],rt)
write('raw_SAT_terminal_checks.csv',['state','max_n','pairs','identical','separate','point','edge','status'],rbase)
print('RAW',len(rstates),len(rt),len(rbase),sum(r[2] for r in rbase),flush=True)
# Four-type atom grammar graph. Base checks examine actual polygons, and classify full atom pairs.
roots,states,arcs,pruned=abuild();sid={st:i for i,st in enumerate(states)};rootrows=[];arcrows=[];leafrows=[];cat=set()
for t,i,j,drop,st in roots:
 a,ii,b,jj,e,r,d=st;out=aoutside(a,ii,b,jj,d);assert out or st in sid
 rootrows.append((t,i,j,drop,str(st),'DISK_PRUNE' if out else sid[st]))
for st in states:
 i=sid[st]
 for ch,drop in achildren(st):
  a,ii,b,jj,e,r,d=ch;out=aoutside(a,ii,b,jj,d);assert out or ch in sid
  arcrows.append((i,drop,str(ch),'DISK_PRUNE' if out else sid[ch]))
 for n in (1,2):
  result=abase(st,n,cat)
  if result is not None:leafrows.append((i,n,sum(result.values()),repr(sorted(result.items())),'PASS_NO_DUPLICATE_NO_PARTIAL_INTERSECTION'))
write('atom_28_states.csv',['id','a','depth_a','b','depth_b','sign','rotation','d'],[(i,*x) for i,x in enumerate(states)])
write('atom_130_sibling_roots.csv',['parent_type','child1','child2','level_drop','state','destination_or_prune'],rootrows)
write('atom_transition_certificate.csv',['parent','level_drop','child','destination_or_prune'],arcrows)
write('atom_SAT_terminal_checks.csv',['state','max_n','contact_atom_pairs','class_counts','status'],leafrows)
# Also directly check assembly bases A/B/K/L with index <=2, and A3/B3 motif catalog.
old={r['signature_SHA256']:r for r in read('prior_10_motifs.csv')};finite=[];counts=Counter();finite_catalog=set()
for n in (1,2,3):
 for t in 'ABKL':
  parts=atoms(t,n);seen=set();contacts=Counter();tile_total=0
  for u,S in parts:
   assert not seen&S;seen|=S;tile_total+=len(S)
  assert seen==actual(t,n)
  for (a,S),(b,T) in combinations(parts,2):
   rel=Counter(classify(x,y) for x in S for y in T)
   assert not any(k for u,k in rel.items() if u=='identical' or u.startswith('BAD'))
   if rel['point'] or rel['edge']:
    k=signature(a,S,b,T);finite_catalog.add(k);h=oldhash(k);assert h in old,(a,b,h)
    contacts[old[h]['type']]+=1
  if n==3 and t in 'AB':
   for h,row in old.items():assert contacts[row['type']]==int(row[t+'3'])
  finite.append((t,n,len(parts),tile_total,sum(contacts.values()),repr(sorted(contacts.items())),'PASS'))
assert len(cat)==len(finite_catalog)==10 and {oldhash(k) for k in cat}==set(old) and cat==finite_catalog
write('motif_hash_match_10.csv',['type','v13_hash','v14_hash','status'],[(old[oldhash(k)]['type'],oldhash(k),oldhash(k),'PASS_EXACT_CONGRUENCE') for k in sorted(cat,key=lambda x:old[oldhash(x)]['type'])])
write('10_full_exact_canonical_motifs.csv',['type','full_two_atom_coordinate_signature'],[(old[oldhash(k)]['type'],repr(k)) for k in sorted(cat,key=lambda x:old[oldhash(x)]['type'])])
write('finite_assembly_bases.csv',['kind','n','atoms','tiles','contact_atom_pairs','motifs','status'],finite)
print('ATOMS',len(states),len(arcrows),len(leafrows),'catalog',len(cat),flush=True)
# Four grammar symbols, independently reconstructed in normalized coordinates.
write('normalized_grammar.csv',['parent','child_id','child_type','depth','sign','rotation','center_coefficient'],[(t,i,*x) for t in 'ABKL' for i,x in enumerate(GR[t])])
write('subset_10_identities.csv',['gap','retained_child','child','source_child','source'],[r for t in 'KL' for r in subset(t)])
write('coverage_22_identities.csv',['gap','source_child','source','case','witness_path','target'],[r for t in 'KL' for r in coverage(t)])
# Finite base of coverage and all exact raw sets up to n=4 (includes a1,a2,b1,b2).
base_cov=[]
for n in range(1,5):
 for t in 'ABKL':
  G=actual(t,n)
  if t in 'AB':R=RAW[t,n]
  elif t=='L':R=RAW['A',n]-trans(RAW['B',n],mul(neg(pp(n)),inv),-1)
  else:R=set() if n==1 else RAW['B',n-1]-trans(RAW['A',n],neg(pp(n-1)),1)
  assert G==R
  base_cov.append((t,n,len(G),len(R),'PASS_EXACT_SET_EQUALITY'))
write('raw_grammar_finite_equivalence.csv',['kind','n','grammar_tiles','raw_or_mask_tiles','status'],base_cov)
# Verify base contacts and containment at n=1,2 independently: proof uses symbolic coverage for n>=3.
# For n=3 all source expansion references have valid indices; only empty K1 may be suppressed.
# n=3 raw mask grandchildren used in coverage may refer to B_{0}? the matched paths all depth <=2.
# Cross-check unquotiented descriptor graph's 285 states.
unquotiented,uarcs,upruned=independent_graph();assert {acanon(*st) for st in unquotiented}==set(states)
write('independent_unquotiented_states.csv',['state','canonical_state_id'],[(str(st),sid[acanon(*st)]) for st in sorted(unquotiented)])
indcat=set();indbase=[]
for ui,st in enumerate(sorted(unquotiented)):
 for n in (1,2):
  res=abase(st,n,indcat)
  if res is not None:indbase.append((ui,sid[acanon(*st)],n,sum(res.values()),repr(sorted(res.items())),'PASS'))
assert indcat==cat
write('independent_unquotiented_517_base_checks.csv',['state_index','canonical_id','max_n','contact_atom_pairs','class_counts','status'],indbase)
print('UNQUOTIENTED',len(unquotiented),uarcs,upruned,len(indbase),flush=True)
# 1000 distinct physical scales. Re-evaluate exact placement arithmetic with independent descriptor composition.
# Select one sibling-root and one graph state per round; all outgoing child cases of that state are recomputed.
# Re-evaluate the radius inequality in absolute coordinates for a rotating sample of pruned arcs.
# No attempts to enumerate astronomically large crystals.
all_pruned=[(st,ch,drop) for st in states for ch,drop in achildren(st) if aoutside(ch[0],ch[1],ch[2],ch[3],ch[-1])]
rounds=[];negative=[];v=mul(s,s) # s^(3-1)
for n in range(3,1003):
 st=states[(n-3)%len(states)];a,i,b,j,e,r,d=st
 X=desc(a,i);Y=desc(b,j,e,r,d)
 L=expand(X) if i==0 else [X];R=expand(Y) if j==0 else [Y]
 # compare *absolute* centers with normalized child coordinates at remaining level.
 cnt=0
 for U0,V0 in product(L,R):
  q=pair(U0,V0);a0,ii,f,k,c=U0;b0,jj,g,l,h=V0;drop=min(ii,jj)
  vv=v
  for _ in range(drop):vv=mul(vv,si)
  absd=rot(sub(mul(h,v),mul(c,v)),f,-k%5)
  assert absd==mul(q[-1],vv)
  assert acanon(*q) in sid or aoutside(q[0],q[1],q[2],q[3],q[-1]);cnt+=1
 parent,ch,drop=all_pruned[(n-3)%len(all_pruned)];a0,i0,b0,j0,e0,r0,d0=ch
 vv=v
 for _ in range(drop):vv=mul(vv,si)
 dist=dot2(mul(d0,vv),mul(d0,vv));rad=qadd(arad(a0,i0),arad(b0,j0));squared=qmul(rad,rad);qscale=dot2(vv,vv)
 threshold=qmul(squared,qscale);assert sign(qsub(dist,threshold))>0
 # perturb exact transformation at this scale (not a random geometric crystal counterexample).
 mutant=plus(absd,I);assert mutant!=mul(q[-1],mul(v,spow_neg(min(U0[1],V0[1]))))
 negative.append((n,'translation_plus_one',str(st),'REJECTED'))
 rounds.append((n,sid[st],cnt,1,1,'PASS'))
 v=mul(v,s)
write('1000_distinct_scale_rounds.csv',['n','state_id','absolute_child_placement_checks','absolute_radius_checks','negative_checks','status'],rounds)
write('1000_mutated_placements.csv',['n','mutation','state','status'],negative)
summary={'raw_states':len(rstates),'raw_arcs':len(rt),'raw_pruned':rpruned,'raw_terminal_cases':len(rbase),'raw_terminal_tile_pairs':sum(r[2] for r in rbase),'atom_states':len(states),'atom_roots':len(roots),'atom_root_pruned':sum(r[-1]=='DISK_PRUNE' for r in rootrows),'atom_arcs':len(arcrows),'atom_pruned':pruned,'atom_terminal_cases':len(leafrows),'unquotiented_states':len(unquotiented),'unquotiented_terminal_cases':len(indbase),'unquotiented_arcs':uarcs,'unquotiented_pruned':upruned,'motif_classes':len(cat),'scale_rounds':len(rounds),'absolute_placement_checks':sum(r[2] for r in rounds),'prior_manifest_matches':len(priorchecks)}
(ROOT/'RESULTS.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8');print(summary,flush=True)
