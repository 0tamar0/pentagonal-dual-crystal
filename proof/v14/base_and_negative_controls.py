"""Independent exhaustive finite base and explicit bad-intersection controls."""
from automaton import *
from itertools import combinations

def write(name,header,rows):
 with (ROOT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.writer(f);w.writerow(header);w.writerows(rows)
rows=[];data=gen(3)
for n in (1,2,3):
 for t in 'AB':
  S=data[t,n];C=Counter(classify(a,b) for a,b in combinations(S,2));assert not any(v for k,v in C.items() if k.startswith('BAD') or k=='identical')
  rows.append((t,n,len(S),sum(C.values()),C['separate'],C['point'],C['edge'],'PASS'))
assert sum(r[3] for r in rows)==105990
write('105990_independent_base_polygon_SAT.csv',['crystal','n','tiles','all_unordered_pairs','separate','point','edge','status'],rows)
# polygon face tests must reject geometrically invalid intersections, not just wrong coefficients.
A=(Z,1);B=(neg(I),-1);assert classify(A,B)=='edge'
# A and B share the support line normal 1; perturb within the line by a short nonzero tangent.
tangent=mul(sub(zz[1],zz[4]),mul(mul(inv,inv),inv));partial=(plus(B[0],tangent),B[1]);assert classify(A,partial)=='BAD_partial_edge'
# vertex versus interior point of edge (golden ratio interpolation along edge, not endpoint).
face=[v for v in verts(A) if dot2(v,I)==interval(A,I)[0]];assert len(face)==2
p=plus(face[0],mul(inv,sub(face[1],face[0])))
victim=(sub(p,inv),1);assert classify(A,victim)=='BAD_partial_point'
area=(Z,-1);assert classify(A,area)=='BAD_area'
write('explicit_geometric_negative_controls.csv',['case','test_polygon1','test_polygon2','classification','status'],[(k,str(A),str(X),classify(A,X),'REJECTED') for k,X in [('partial_edge',partial),('vertex_inside_edge',victim),('positive_area',area)]])
print('105990 exact polygon pairs, 3 explicit invalid polygon contacts rejected',flush=True)
