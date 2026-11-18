"""Independent exact finite-state certificate for pentagon intersection safety.
No previously proved mask or contact theorem is used. Python standard library.
"""
from collections import deque,Counter
from functools import lru_cache
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parent
Z=(0,0,0,0); I=(1,0,0,0); z=(0,1,0,0)
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def times(a,k):return tuple(x*k for x in a)
def sub(a,b):return plus(a,neg(b))
def rz(a):
 a,b,c,d=a;return -d,a-d,b-d,c-d
def mul(a,b):
 r=Z
 for k in b:r=plus(r,times(a,k));a=rz(a)
 return r
zz=[I,z,rz(z),rz(rz(z)),rz(rz(rz(z)))]; phi=neg(plus(zz[2],zz[3])); ph2=plus(I,phi); inv=sub(phi,I); s=neg(ph2)
def conj(a):
 x,y,z,w=a;return x-y,-y,w-y,z-y
# real part doubled, coefficients of 1 and phi
def re2(a):
 x,y,z,w=a;return 2*x-y,y-z-w
def sign(q):
 a,b=q; t=2*a+b
 if not b:return (t>0)-(t<0)
 if b>0:return 1 if t>=0 or 5*b*b>t*t else -1
 return -1 if t<=0 or 5*b*b>t*t else 1
def dot2(a,b):return re2(mul(a,conj(b)))
def qsub(a,b):return a[0]-b[0],a[1]-b[1]
def qadd(a,b):return a[0]+b[0],a[1]+b[1]
def qmul(a,b):return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
# radius A normalized phi, B normalized phi^2. For a level one smaller multiply phi^-2.
RAD={('A',0):(0,1),('B',0):(1,1),('A',1):(-1,1),('B',1):(1,0)}
@lru_cache(None)
def canon(a,i,b,j,e,d):
 # sort by depth then type; first side sign +. Reflection and 5-fold rotational invariance.
 if (i,a)>(j,b):a,i,b,j,d=b,j,a,i,times(d,-e)
 candidates=[]
 for t in (d,conj(d)):
  for _ in range(5):candidates.append(t);t=rz(t)
 return a,i,b,j,e,min(candidates)
@lru_cache(None)
def outside(a,i,b,j,d):
 r=qadd(RAD[a,i],RAD[b,j]); lim=qmul(r,r); lim=(lim[0]*2,lim[1]*2)
 return sign(qsub(dot2(d,d),lim))>0
# branch data = (type, depth relative to m, orientation relative to parent, translation / s^(m-1))
def branches(a,i):
 if i:return [(a,1,1,Z)]
 rad=I if a=='A' else phi
 out=[(a,1,1,Z)]+[(a,1,-1,neg(mul(rad,u))) for u in zz]
 if a=='A':out += [('B',2,-1,neg(mul(inv,u))) for u in zz]
 else:out += [('A',1,1,neg(u)) for u in zz]
 return out
BR={(a,i):branches(a,i) for a in 'AB' for i in (0,1)}
def children(st):
 a,i,b,j,e,d=st
 for a1,i1,f,c in BR[a,i]:
  for b1,j1,g,h in BR[b,j]:
   # when both children depth=2, max drops 2, relative offset multiplies s^2.
   drop=min(i1,j1); assert drop in (1,2)
   D=times(plus(sub(d,c),times(h,e)),f)
   for _ in range(drop):D=mul(D,s)
   yield canon(a1,i1-drop,b1,j1-drop,e*f*g,D),drop

def build():
 roots=[canon(a,0,a,0,1,Z) for a in 'AB']
 seen=set(roots); queue=deque(roots); arcs=0; pruned=0
 while queue:
  st=queue.popleft()
  for ch,drop in children(st):
   arcs+=1
   a,i,b,j,e,d=ch
   if outside(a,i,b,j,d):pruned+=1;continue
   if ch not in seen:seen.add(ch);queue.append(ch)
  if len(seen)>200000:raise RuntimeError('state limit')
 return roots,sorted(seen),arcs,pruned
if __name__=='__main__':
 roots,states,arcs,pruned=build();print('states',len(states),'arcs',arcs,'pruned',pruned,flush=True)
 print(Counter((a,i,b,j) for a,i,b,j,e,d in states))

def trans(S,c,e):return {(plus(c,times(p,e)),e*f) for p,f in S}
def gen(n):
 A={(Z,1)}|{(neg(u),-1) for u in zz};B={(neg(u),1) for u in zz};allsets={('A',1):A,('B',1):B}
 v=I
 for k in range(2,n+1):
  v=mul(v,s);new={}
  for x in 'AB':
   T=set()
   for y,j,e,c in BR[x,0]:
    if k-j<1:continue
    T|=trans(allsets[y,k-j],mul(c,v),e)
   new[x,k]=T
  allsets.update(new)
 return allsets
@lru_cache(None)
def verts(tile):
 p,e=tile;return tuple(plus(p,times(mul(inv,u),e)) for u in zz)
@lru_cache(None)
def interval(tile,u):
 vals=[dot2(v,u) for v in verts(tile)]
 from functools import cmp_to_key
 key=cmp_to_key(lambda a,b:sign(qsub(a,b)))
 return min(vals,key=key),max(vals,key=key)
def qmin(a,b):return a if sign(qsub(a,b))<=0 else b
def qmax(a,b):return a if sign(qsub(a,b))>=0 else b
def classify(a,b):
 if a==b:return 'identical'
 # every edge-normal is among +/- zeta^j, as both orientations are +/-.
 tangent_axis=None
 for u in zz:
  alo,ahi=interval(a,u);blo,bhi=interval(b,u)
  cmp=sign(qsub(qmin(ahi,bhi),qmax(alo,blo)))
  if cmp<0:return 'separate'
  if cmp==0:tangent_axis=(u,ahi if ahi==blo else bhi)
 if tangent_axis is None:return 'BAD_area'
 u,h=tangent_axis
 vv=[v for v in verts(a) if dot2(v,u)==h];ww=[v for v in verts(b) if dot2(v,u)==h]
 assert vv and ww
 t=mul(u,sub(zz[1],zz[4]))
 va=[dot2(v,t) for v in vv];vb=[dot2(v,t) for v in ww]
 al=va[0];ah=va[0];bl=vb[0];bh=vb[0]
 for v in va:al=qmin(al,v);ah=qmax(ah,v)
 for v in vb:bl=qmin(bl,v);bh=qmax(bh,v)
 lo=qmax(al,bl);hi=qmin(ah,bh);sg=sign(qsub(hi,lo))
 if sg<0:return 'separate' # redundant if SAT accurate
 if sg==0:return 'point' if lo in va and lo in vb else 'BAD_partial_point'
 return 'edge' if al==bl==lo and ah==bh==hi and len(va)==len(vb)==2 else 'BAD_partial_edge'
def base(st,n,sets):
 a,i,b,j,e,d=st
 if n-max(i,j)<1:return None
 scale=I if n==1 else s
 X=sets[a,n-i];Y=trans(sets[b,n-j],mul(d,scale),e)
 return Counter(classify(x,y) for x in X for y in Y)
if __name__=='__main__':
 sets=gen(2);print({k:len(v) for k,v in sets.items()})
 for i,st in enumerate(states):
  print(i,st,[base(st,n,sets) for n in (1,2)],flush=True)
