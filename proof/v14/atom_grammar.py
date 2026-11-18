"""Four-type grammar: all translations are normalized by s^(n-1).
K = a in the B-insertion coordinates, L = b in A-insertion coordinates.
"""
from automaton import *
from itertools import combinations
# transform using signed rotation pair: j in Z/5, e in {-1,1}.
def rot(p,e,j):return times(mul(p,zz[j]),e)
def mapset(S,c,e,j):return {(plus(c,rot(p,e,j)),f*e) for p,f in S}
si=neg(mul(inv,inv)) # s^-1
assert mul(si,s)==I
# child records (type,depth,orientation,rotation,translation)
def grammar(t):
 if t in 'AB':
  rad=I if t=='A' else phi
  out=[(t,1,1,0,Z)]+[(t,1,-1,0,neg(mul(rad,u))) for u in zz]
  if t=='A':out += [('K',1,-1,j,neg(mul(inv,u))) for j,u in enumerate(zz)]
  else:out += [('L',1,1,j,neg(u)) for j,u in enumerate(zz)]
 elif t=='L':
  out=[('A',1,1,0,Z),('A',1,-1,0,neg(I))]+[('L',1,-1,0,neg(zz[j])) for j in (1,4)]+[('K',1,-1,0,neg(inv))]
 else:
  # B_(n-1) child centers are v_(n-1) * phi = -s^(n-2)*phi
  w=neg(mul(si,phi));v=neg(si)
  out=[('B',2,1,0,Z),('B',2,-1,0,w)]+[('K',1,-1,0,mul(w,zz[j])) for j in (1,4)]+[('L',2,1,0,v)]
 return out
GR={t:grammar(t) for t in 'ABKL'}
RAW=gen(4)
# actual gap sets from *raw* at bases
K2=RAW['B',1]-trans(RAW['A',2],neg(s),1) # v2=-s
L1=RAW['A',1]-trans(RAW['B',1],mul(neg(s),inv),-1)
assert len(K2)==1 and len(L1)==4
BASE={('A',1):RAW['A',1],('B',1):RAW['B',1],('K',1):set(),('K',2):K2,('L',1):L1}

def children_for(t,n):
 if (t,n) in BASE:return []
 return [x for x in GR[t] if n-x[1]>=1 and not(x[0]=='K' and n-x[1]==1)]
@lru_cache(None)
def pp(n):
 x=I
 for _ in range(n):x=mul(x,s)
 return x
@lru_cache(None)
def actual(t,n):
 if (t,n) in BASE:return BASE[t,n]
 parts=[mapset(actual(y,n-d),mul(c,pp(n-1)),e,j) for y,d,e,j,c in children_for(t,n)]
 assert sum(map(len,parts))==len(set().union(*parts)),(t,n,'duplicate grammar')
 return set().union(*parts)
if __name__=='__main__':
 for n in range(1,5):
  for t in 'ABKL':
   a=actual(t,n)
   if t in 'AB':b=RAW[t,n]
   elif t=='L':b=RAW['A',n]-trans(RAW['B',n],mul(neg(pp(n)),inv),-1)
   else:b=set() if n==1 else RAW['B',n-1]-trans(RAW['A',n],neg(pp(n-1)),1)
   print(t,n,len(a),len(b),a==b,flush=True)
   assert a==b
@lru_cache(None)
def acanon(a,i,b,j,e,r,d):
 if (i,a)>(j,b):
  a,i,b,j,d,r=b,j,a,i,rot(d,-e,-r%5),-r%5
 opts=[]
 for k in range(5) if a in 'AB' else (0,):
  D=rot(d,1,k);R=(r+k)%5 if b not in 'AB' else 0
  opts.append((R,D));opts.append((-R%5,conj(D)))
 r,d=min(opts)
 return a,i,b,j,e,r,d

def arad(t,i):
 r={'A':(0,1),'B':(1,1),'K':(1,0),'L':(0,1)}[t]
 for _ in range(i):r=qmul(r,(2,-1)) # phi^-2
 return r
@lru_cache(None)
def aoutside(a,i,b,j,d):
 r=qadd(arad(a,i),arad(b,j));lim=qmul(r,r)
 return sign(qsub(dot2(d,d),times(lim,2)))>0

def achildren(st):
 a,i,b,j,e,r,d=st
 for aa,ii,f,k,c in (GR[a] if i==0 else [(a,i,1,0,Z)]):
  for bb,jj,g,l,h in (GR[b] if j==0 else [(b,j,1,0,Z)]):
   drop=min(ii,jj);assert drop in (1,2)
   D=rot(plus(sub(d,c),rot(h,e,r)),f,-k%5)
   for _ in range(drop):D=mul(D,s)
   yield acanon(aa,ii-drop,bb,jj-drop,e*f*g,(r+l-k)%5,D),drop

def aroots():
 out=[]
 for t in 'ABKL':
  for i,j in combinations(range(len(GR[t])),2):
   a,ii,f,k,c=GR[t][i];b,jj,g,l,h=GR[t][j];drop=min(ii,jj)
   D=rot(sub(h,c),f,-k%5)
   for _ in range(drop):D=mul(D,s)
   st=acanon(a,ii-drop,b,jj-drop,f*g,(l-k)%5,D)
   out.append((t,i,j,drop,st))
 return out

def abuild():
 roots=aroots();seen=set();queue=deque();pruned=0;arcs=0
 for *_,st in roots:
  a,i,b,j,e,r,d=st
  if not aoutside(a,i,b,j,d) and st not in seen:seen.add(st);queue.append(st)
 while queue:
  st=queue.popleft()
  for ch,drop in achildren(st):
   arcs+=1;a,i,b,j,e,r,d=ch
   if aoutside(a,i,b,j,d):pruned+=1;continue
   if ch not in seen:seen.add(ch);queue.append(ch)
  if len(seen)>100000:raise RuntimeError('state limit')
 return roots,sorted(seen),arcs,pruned

@lru_cache(None)
def atoms(t,n):
 if (t,n) in BASE:return [(t,BASE[t,n])] if BASE[t,n] else []
 out=[]
 for y,d,e,j,c in children_for(t,n):
  for a,S in atoms(y,n-d):out.append((a,mapset(S,mul(c,pp(n-1)),e,j)))
 return out

def signature(a,X,b,Y):
 if a>b:a,b,X,Y=b,a,Y,X
 opts=[]
 for mirror in (0,1):
  for e in (1,-1):
   for j in range(5):
    def F(S):return sorted((rot(conj(p) if mirror else p,e,j),f*e) for p,f in S)
    M,N=F(X),F(Y);p=min(q for q,f in M+N)
    U,V=tuple((sub(q,p),f) for q,f in M),tuple((sub(q,p),f) for q,f in N)
    if a==b and U>V:U,V=V,U
    opts.append((a,b,U,V))
 return min(opts)

def abase(st,n,catalog):
 a,i,b,j,e,r,d=st
 if n-max(i,j)<1:return None
 X=atoms(a,n-i);Y=[(u,mapset(S,mul(d,pp(n-1)),e,r)) for u,S in atoms(b,n-j)]
 cnt=Counter()
 for u,S in X:
  for v,T in Y:
   rel=Counter(classify(x,y) for x in S for y in T)
   bad=sum(k for t,k in rel.items() if t=='identical' or t.startswith('BAD'))
   if bad:print('BAD',st,n,u,v,rel,flush=True);raise AssertionError('bad atom contact')
   if rel['edge'] or rel['point']:
    sig=signature(u,S,v,T)
    catalog.add(sig);cnt[(u,v,'edge' if rel['edge'] else 'point')]+=1
 return cnt
if __name__=='__main__':
 rr,ss,ar,pr=abuild();print('ATOM states',len(ss),'arcs',ar,'pruned',pr,'roots',len(rr),flush=True)
 cat=set()
 for i,st in enumerate(ss):
  for n in (1,2):abase(st,n,cat)
 print('ATOM catalog from state bases',len(cat),Counter((a,b) for a,b,_,_ in cat),flush=True)
 # Include sibling contacts at n=2 because A2,B2,L2 may generate no K1 and n=2 root normalization can have max1.
 for n in (1,2,3):
  for t in 'ABKL':
   X=atoms(t,n)
   for (a,S),(b,T) in combinations(X,2):
    rel=Counter(classify(x,y) for x in S for y in T)
    assert not any(k for u,k in rel.items() if u=='identical' or u.startswith('BAD'))
    if rel['point'] or rel['edge']:cat.add(signature(a,S,b,T))
 print('ATOM catalog with finite assembly bases',len(cat),Counter((a,b) for a,b,_,_ in cat),flush=True)
