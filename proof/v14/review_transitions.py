"""Independent reconstruction of all previously proposed standard contact rules.
No old geometric lemmas or proof source are imported. Uses only v14 graph.
"""
from independent import *

def write(name,header,rows):
 with (ROOT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.writer(f);w.writerow(header);w.writerows(rows)
roots,states,_,_=abuild();sid={st:i for i,st in enumerate(states)}
active={st for st in states if any(abase(st,n,set()) for n in (1,2))}
while True:
 new=active|{st for st in states if any(ch in active for ch,d in achildren(st))}
 if new==active:break
 active=new
assert len(active)==13
# No active terminal can be reached from an inactive vertex; this is an inductive certificate.
for st in set(states)-active:
 assert not any(ch in active for ch,d in achildren(st))
 assert all(not abase(st,n,set()) for n in (1,2))
write('28_states_contact_activity.csv',['id','state','status'],[(sid[st],str(st),'ACTIVE' if st in active else 'ALL_GENERATIONS_NO_CONTACT') for st in states])
# Direct normalized standard geometry, not selected from numbers.
stdraw={
'AA_edge':('A',0,'A',0,-1,0,ph2),
'AA_point':('A',0,'A',0,1,0,mul(ph2,sub(I,z))),
'H':('A',0,'B',0,1,0,s),
'P':('B',0,'A',1,-1,0,phi),
'D':('A',0,'B',1,-1,0,phi),
'C':('A',0,'K',0,-1,0,phi),
'O':('A',0,'K',0,1,0,I),
'Bb_c':('B',0,'L',0,1,0,ph2),
'Bb_o':('B',0,'L',0,-1,0,phi),
'Ab_edge':('A',0,'L',0,-1,0,mul(ph2,z)),
'Ab_point':('A',0,'L',0,1,0,sub(ph2,mul(ph2,z))),
'bb':('L',0,'L',0,1,1,mul(ph2,sub(z,I))),
}
std={k:acanon(*st) for k,st in stdraw.items()}
assert {sid[st] for k,st in std.items() if k!='P'}=={0,1,3,4,6,8,9,10,18,19,26}
name={st:k for k,st in std.items()};name[states[11]]='Ab_diff_point';name[states[12]]='Ab_diff_edge'
assert set(name)==active|{std['P']}
expected={'AA_edge':{'AA_edge':2},'AA_point':{'AA_point':1},'H':{'H':4,'Ab_edge':2},'P':{'H':2,'Ab_edge':2},'D':{'H':2},'C':{'C':2,'D':2,'Ab_diff_edge':2},'O':{'O':2,'D':2,'Ab_diff_point':2},'Bb_c':{'Bb_c':2,'H':2,'Ab_edge':2},'Bb_o':{'Bb_o':2,'H':2,'Ab_edge':2},'Ab_edge':{'Ab_edge':1},'Ab_point':{'Ab_point':1},'bb':{'bb':1},'Ab_diff_edge':{'Ab_edge':1},'Ab_diff_point':{'Ab_point':1}}
transition=[];counts={}
for st,k in name.items():
 # This also checks the extra P root, which is not needed for the 28-state closure of assembly siblings.
 children_=list(achildren(st));assert all(ch in sid or aoutside(ch[0],ch[1],ch[2],ch[3],ch[-1]) for ch,d in children_)
 c=Counter(name[ch] for ch,d in children_ if ch in active)
 assert all(d==1 for ch,d in children_ if ch in active)
 assert c==expected[k],(k,c,expected[k]);counts[k]=c
 transition.append((k,sid.get(st,'EXTERNAL_P'),str(st),repr(sorted(c.items())),len(children_),sum(c.values()),'ALL_GENERATIONS_N>=3; BASES_SEPARATE'))
write('14_standard_and_auxiliary_rules.csv',['rule','id','state','nonzero_child_rules','all_children','active_children','status'],sorted(transition))
# Group B-b transitions in precisely the original 2+1 form: hold B fixed while splitting L.
group=[]
for k in ('Bb_c','Bb_o'):
 st=std[k];a,i,b,j,e,r,d=st;X=desc(a,i);Y=desc(b,j,e,r,d)
 P=0;lower=0
 for q,U in enumerate(expand(Y)):
  p=acanon(*pair(X,U));children_=list(achildren(p));assert all(ch in sid or aoutside(ch[0],ch[1],ch[2],ch[3],ch[-1]) for ch,g in children_)
  c=Counter(name[ch] for ch,g in children_ if ch in active)
  assert p[1]==0 and p[3]==1
  if p==std['P']:P+=1;assert c==counts['P'];tag='P'
  elif c==Counter({k:1}):lower+=1;tag=k+' previous'
  else:assert not c;tag='NONE'
  # max generation=2 is the only extra base; K1 empty, all other groups checked directly.
  base_rel=abase(p,2,set());assert bool(base_rel)==(tag!='NONE')
  group.append((k,q,str(p),tag,repr(sorted(c.items())),repr(sorted(base_rel.items())),'PASS'))
 assert (P,lower)==(1,2)
write('Bb_2_plus_1_exact_grouping.csv',['rule','gap_child','pair_state','branch','expanded_active_children','base_n2','status'],group)
# Active sibling contacts inside each assembly, verified before the all-gen proof.
siblings=[]
for t,i,j,drop,st in roots:
 if st in active:siblings.append((t,i,j,drop,name[st],str(st)))
write('active_sibling_root_rules.csv',['parent','child1','child2','drop','standard_rule','state'],siblings)
# Direct tile-edge incidence base counts (using exact SAT; allowed edge intersection = entire common edge).
def ecount(st,n):
 a,i,b,j,e,r,d=st
 if n<max(i,j)+1:return 0
 X=actual(a,n-i);Y=mapset(actual(b,n-j),mul(d,pp(n-1)),e,r)
 return sum(classify(x,y)=='edge' for x in X for y in Y)
V1={k:ecount(st,1) for st,k in name.items()};V2={k:ecount(st,2) for st,k in name.items()}
# n=2 transition is valid for A/B/L; K2 is atomic and C/O at n2 are direct bases.
# Check exact first step for all other nonzero-domain states (incl. P/D).
for k in expected:
 if k in ('C','O'):continue
 assert V2[k]==sum(m*V1[l] for l,m in counts[k].items()),(k,V2[k],counts[k],V1)
# Check that each terminal base actually yields the required atom signatures (tags in v14 catalog).
rows=[]
for n in range(3,1003):
 V={k:sum(m*V2[l] for l,m in c.items()) for k,c in counts.items()}
 h=2*(4**n-1)//3;prev=2*(4**(n-1)-1)//3
 assert V['H']==h and V['Bb_c']==V['Bb_o']==h
 assert V['AA_edge']==2**n and V['AA_point']==0 and V['bb']==0
 assert V['D']==2*prev and V['P']==2*prev+2
 assert V['C']==(4**n+2)//3-2**n and V['O']==(4**n+8)//3-2**(n+1)
 assert V['Ab_edge']==V['Ab_diff_edge']==1 and V['Ab_point']==V['Ab_diff_point']==0
 rows.append((n,V['AA_edge'],V['H'],V['P'],V['D'],V['C'],V['O'],V['Bb_c'],V['Bb_o'],'PASS'))
 V2=V
write('1000_standard_rule_edge_count_regressions.csv',['max_n','AA','H','P_index_n_minus_1','D_index_n_minus_1','C','O','Bb_c','Bb_o','status'],rows)
print('14 rules, 2+1 grouping, active siblings',len(siblings),'and 1000 DP index checks PASS',flush=True)
