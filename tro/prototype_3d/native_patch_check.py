"""Regression checks for native covers and radius-dependent derivative certificates."""
import os
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from scipy.optimize import brentq
from proto3d import Spline3,surface_boxes
import summed as S
import summed3d as S3

f=Spline3(np.full(3,-.08),np.full(3,.08),.02).fit(lambda p:np.linalg.norm(p,axis=1)-.045,.01)
lo,side=surface_boxes(f)
assert side==f.h
assert np.allclose((lo-f.lo)/f.h,np.round((lo-f.lo)/f.h))
body=S.prepare_body(f,0.,lo+side/2,side,3)
rng=np.random.default_rng(53)
directions=rng.normal(size=(80,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
for d in directions:
    root=brentq(lambda t:float(f.eval((t*d)[None])[0]),.025,.065)
    p=root*d
    assert np.any(np.all((p>=lo-1e-12)&(p<=lo+side+1e-12),axis=1))
cell=S.third_cellwise(f)
B={'f':f}
S3.prep_field(B,.003)
S3.prep_field(B,body['r'])
assert B['_M3_reach']>=body['r']
points=rng.uniform(-.03,.03,(100,3))
exact_cell_max=S.box_bound(f,cell,points,body['r'])
assert np.all(B['M3'].eval(points)>=exact_cell_max-1e-10)
large_reach=1.6*f.h
S3.prep_field(B,large_reach)
hessian_max=S.box_bound(f,f.hessian_bound(cellwise=True),points,large_reach)
assert B['_M2_reach']>=large_reach
assert np.all(B['M2_radius'].eval(points)>=hessian_max-1e-10)
assert S.REFINE_DEPTH==0
print('PASS: native cell alignment, 80 surface roots covered, full-patch Taylor preparation,')
print('      enlarged cached radius bounds every intersected cell at 100 centers; refinement default is zero.')
print('      Hessian pruning bounds also cover patches wider than the other field cells.')
