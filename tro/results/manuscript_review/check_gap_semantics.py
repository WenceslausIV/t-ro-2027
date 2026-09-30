"""Regression: a negative Taylor distance bound is not a collision witness."""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'prototype_3d'))
import franka3d as F
from refine_collision_audit import audit_pair


def main():
    obstacle = dict(boxes=[((0.,0.,0.),(.05,.05,.05))],
                    blo=np.full(3,-.05), bhi=np.full(3,.05))
    original = F.fk
    F.fk = lambda q: ([np.eye(4)], None, None)
    results = {}
    try:
        for name, height in [('separated',.051),('intersecting',.049),('far',.3)]:
            vertices=np.array([[-.001,-.001,height],[.001,-.001,height],[0,.002,height]])
            center=vertices.mean(axis=0)
            link=dict(frame=0,V=vertices,gt=center[None],
                      gt_r=np.array([np.linalg.norm(vertices-center,axis=1).max()]))
            lower,witness=F.real_gap(np.zeros(7),[link],{'box':obstacle},return_witness=True)
            adaptive=audit_pair(vertices[None],obstacle)
            results[name]=dict(lower=lower,witness=witness,adaptive=adaptive)
        assert results['separated']['lower'] < 0 < results['separated']['witness']
        assert results['separated']['adaptive']['status']=='separated'
        assert results['intersecting']['witness']<0
        assert results['intersecting']['adaptive']['status']=='collision'
        assert 0<results['far']['lower']<=.25 and np.isfinite(results['far']['lower'])
    finally:
        F.fk=original
    vertices=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    faces=np.array([[0,1,2]])
    points=np.array([[.2,.2,.01],[.2,.2,10.]])
    distances=F.mesh_distance(vertices,faces,points,.1)
    assert np.isclose(distances[0],.01)
    assert 0<distances[1]<=10 and np.isfinite(distances[1])
    assert F.mesh_distance(vertices,faces,points[1:],.1).shape==(1,)
    assert F.mesh_distance(vertices,faces,np.empty((0,3)),.1).shape==(0,)
    results['mesh_distance_near_far']=distances.tolist()
    Path(__file__).with_suffix('.json').write_text(json.dumps(results,indent=2))
    print('Gap semantics and adaptive triangle checks passed.')


if __name__=='__main__':
    main()
