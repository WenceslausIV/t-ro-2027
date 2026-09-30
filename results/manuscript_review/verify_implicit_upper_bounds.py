"""Verify a saved boundary level using only exhaustive upper-bound/exclusion tests.

No projected lower bound enters this audit. Refinement is chunked to limit temporary storage.
An exhausted work budget is reported as unresolved, never as a failed enclosure.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
sys.path.insert(0, str(ROOT / 'prototype_3d/interactive'))
import franka3d as F
from proto3d import sdf_boxes
import server as S

OUT = Path(__file__).resolve().parent


def verify(name, body, sdf, max_depth=12, max_boxes=4000000):
    f = F.spline_from(body['f']); M = F.hessian_majorant(f)
    side = min(.008, 2*f.h/np.sqrt(3))
    axes = [np.arange(body['blo'][a]-side/2, body['bhi'][a]+side, side) for a in range(3)]
    centers = np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3)
    radius = side*np.sqrt(3)/2
    centers = centers[np.abs(sdf(centers)) <= radius]
    history, minimum_margin = [], np.inf
    for depth in range(max_depth+1):
        if not len(centers):
            return dict(status='verified', level=body['l'],minimum_discard_margin=minimum_margin,history=history)
        pending = []
        for start in range(0,len(centers),50000):
            C = centers[start:start+50000]
            v, g = f.eval(C,order=1)
            upper = v+np.linalg.norm(g,axis=1)*radius+.5*M.eval(C)*radius**2
            keep = upper >= body['l']
            if (~keep).any():
                minimum_margin = min(minimum_margin,float((body['l']-upper[~keep]).min()))
            pending.append(C[keep])
        C = np.concatenate(pending)
        history.append(dict(depth=depth,side=side,evaluated=len(centers),remaining=len(C)))
        print(name,history[-1],flush=True)
        if not len(C):
            return dict(status='verified', level=body['l'],minimum_discard_margin=minimum_margin,history=history)
        if depth == max_depth:
            break
        side /= 2; radius = side*np.sqrt(3)/2
        children = []
        for start in range(0,len(C),50000):
            block = (C[start:start+50000,None,:] + (F.CORNERS*side)[None]).reshape(-1,3)
            children.append(block[np.abs(sdf(block)) <= radius])
            if sum(len(p) for p in children) > max_boxes:
                return dict(status='unresolved_box_limit',level=body['l'],history=history)
        centers = np.concatenate(children)
    return dict(status='unresolved_depth_limit',level=body['l'],history=history)


def main():
    old = np.load(F.CACHE,allow_pickle=True)['d'].item()
    objects = np.load(ROOT/'prototype_3d/interactive/cache_objects.npz',allow_pickle=True)['d'].item()
    requests = [('arch',old['obst']['arch'],lambda p:sdf_boxes(p,F.OBSTACLES['arch'])),
                ('rounded_star',objects['star tube (rounded)'],S.shape_spec('star tube (rounded)')['sdf'])]
    destination = OUT/'verify_implicit_upper_bounds.json'
    result = json.loads(destination.read_text()) if destination.exists() else {}
    for name, body, sdf in requests:
        if result.get(name, {}).get('status') == 'verified':
            continue
        result[name] = verify(name,body,sdf)
        destination.write_text(json.dumps(result,indent=2))
    assert all(r['status']=='verified' for r in result.values()), result


if __name__=='__main__':
    main()
