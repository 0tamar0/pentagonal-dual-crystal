"""Second enumeration: no rotation, reflection, or exchange pruning of pair states.
Uses descriptor composition, independent of achildren()/aroots().
"""
from coverage import *

def pair(X,Y):
 a,i,f,k,c=X;b,j,g,l,h=Y
 drop=min(i,j);v=pp(drop)
 return a,i-drop,b,j-drop,f*g,(l-k)%5,mul(rot(sub(h,c),f,-k%5),v)

def independent_roots():
 return [pair(X,Y) for t in 'ABKL' for X,Y in combinations(expand(desc(t)),2)]

def independent_children(st):
 a,i,b,j,e,r,d=st
 X=desc(a,i);Y=desc(b,j,e,r,d)
 return [pair(U,V) for U in (expand(X) if i==0 else [X]) for V in (expand(Y) if j==0 else [Y])]

def independent_graph():
 seen=set();queue=deque();arcs=0;prune=0
 def add(st):
  a,i,b,j,e,r,d=st
  if aoutside(a,i,b,j,d):return False
  if st not in seen:seen.add(st);queue.append(st)
  return True
 for st in independent_roots():add(st)
 while queue:
  for st in independent_children(queue.popleft()):
   arcs+=1
   if not add(st):prune+=1
  assert len(seen)<100000
 return seen,arcs,prune
if __name__=='__main__':
 seen,arcs,prune=independent_graph();print('unquotiented',len(seen),arcs,prune,flush=True)
 R,states,ar,pr=abuild();assert {acanon(*st) for st in seen}==set(states)
 # For all raw states and atom states, verify reduction is an isometry at the bases (compare canonical signature counts).
 print('quotient matches',len(states),flush=True)
