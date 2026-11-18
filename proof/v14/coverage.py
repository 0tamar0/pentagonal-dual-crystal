"""Exact homogeneous inclusion certificates for the four-type grammar.
Every h/rotation is represented in Z[zeta5] (normalized by s^(n-1)).
"""
from atom_grammar import *
from itertools import product

def desc(t,d=0,e=1,j=0,c=Z):return t,d,e,j,c
@lru_cache(None)
def spow_neg(d):
 out=I
 for _ in range(d):out=mul(out,si)
 return out

def expand(X,raw=False):
 t,d,e,j,c=X
 if raw and t in 'AB':spec=[(u,k,f,0,p) for u,k,f,p in BR[t,0]]
 else:spec=GR[t]
 return [(u,d+k,e*f,(j+r)%5,plus(c,rot(mul(p,spow_neg(d)),e,j))) for u,k,f,r,p in spec]

def sym_eq(X,Y):
 t,d,e,j,c=X;u,k,f,r,p=Y
 return (t,d,e,c)==(u,k,f,p) and (t in 'AB' or j==r)

def envelope(X):
 t,d,e,j,c=X
 if t=='K':return 'B',d+1,e,j,c
 if t=='L':return 'A',d,e,j,c
 return X

def direct(X,Y):return sym_eq(X,Y) or sym_eq(envelope(X),Y)

def depth_list(root,levels,raw=False):
 out=[((),root)];q=[((),root)]
 for _ in range(levels):
  nq=[]
  for path,X in q:
   if raw and X[0] not in 'AB':continue
   for i,Y in enumerate(expand(X,raw)):nq.append((path+(i,),Y))
  out+=nq;q=nq
 return out

def coverage(t):
 # K_n source=B_(n-1), mask=v_n+A_n. L_n source=A_n, mask=h_q(B_n).
 src=desc('B',1) if t=='K' else desc('A')
 mask=desc('A',0,1,0,neg(I)) if t=='K' else desc('B',0,-1,0,phi)
 retained=list(enumerate(expand(desc(t))))
 mm=depth_list(mask,2,True)
 rows=[]
 for i,X in enumerate(expand(src)):
  options=[('RETAIN',str(j),Y) for j,Y in retained if direct(X,Y)]
  options += [('MASK',str(path),Y) for path,Y in mm if direct(X,Y)]
  if options:
   z=options[0];rows.append((t,i,str(X),z[0],z[1],str(z[2])));continue
  # inductive lower coverage source: X = transformed B_(k-1) or A_k.
  for j,Y in retained:
   u,d,e,r,c=Y
   if u not in 'KL' or d<1:continue
   if not sym_eq(X,envelope(Y)):continue
   # transform appropriate complementary mask into target local coordinates
   if u=='K':P=desc('A',d,e,r,plus(c,rot(mul(neg(I),spow_neg(d)),e,r)))
   else:P=desc('B',d,-e,r,plus(c,rot(mul(phi,spow_neg(d)),e,r)))
   for path,Q in mm:
    if sym_eq(P,Q):options.append(('INDUCTIVE_'+u,str(j)+' / '+str(path),Q))
  if not options:raise AssertionError((t,i,X,'no cover'))
  z=options[0];rows.append((t,i,str(X),z[0],z[1],str(z[2])))
 return rows

def subset(t):
 src=desc('B',1) if t=='K' else desc('A')
 yy=expand(src);out=[]
 for i,X in enumerate(expand(desc(t))):
  opts=[(j,Y) for j,Y in enumerate(yy) if direct(X,Y)]
  assert opts,(t,i,X,yy)
  out.append((t,i,str(X),opts[0][0],str(opts[0][1])))
 return out
if __name__=='__main__':
 for t in 'KL':
  print('subset',t,subset(t))
  print('coverage',t,coverage(t))
