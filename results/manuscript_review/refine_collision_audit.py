"""Adaptive all-triangle audit of saved states with inconclusive coarse lower bounds.

An actual mesh point with negative primitive-union defining value witnesses intersection.
A positive centroid-minus-radius bound excludes the whole triangle; otherwise subdivide.
AABB exclusion is geometric. No fitted field is queried. No between-step claim is made.
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prototype_3d'))
import franka3d as F
from proto3d import sdf_boxes

OUT = Path(__file__).resolve().parent


def audit_pair(triangles, obstacle, tol=1e-6):
    boxes, low, high = obstacle['boxes'], obstacle['blo'], obstacle['bhi']
    tmin, tmax = triangles.min(axis=1), triangles.max(axis=1)
    delta = np.maximum(np.maximum(low-tmax, tmin-high), 0)
    aabb = np.linalg.norm(delta, axis=1)
    lower = float(aabb[aabb > 0].min()) if np.any(aabb > 0) else np.inf
    triangles = triangles[aabb == 0]
    for depth in range(15):
        if not len(triangles):
            return dict(status='separated', lower_bound=lower)
        center = triangles.mean(axis=1)
        radius = np.linalg.norm(triangles-center[:, None], axis=2).max(axis=1)
        value = sdf_boxes(center, boxes)
        k = int(np.argmin(value))
        if value[k] < -1e-12:
            return dict(status='collision', witness=center[k].tolist(), witness_value=float(value[k]), depth=depth)
        bounds = value-radius
        safe = bounds > 0
        if safe.any():
            lower = min(lower, float(bounds[safe].min()))
        if safe.all():
            return dict(status='separated', lower_bound=lower)
        if radius[~safe].max() < tol or depth == 14:
            return dict(status='unresolved', lower_bound=float(bounds[~safe].min()), depth=depth)
        a, b, c = (triangles[~safe, j] for j in range(3))
        ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
        triangles = np.concatenate([np.stack(t, axis=1) for t in
                                    ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
    raise AssertionError('unreachable')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*')
    args = ap.parse_args()
    old = np.load(F.CACHE, allow_pickle=True)['d'].item()
    objects = {name: dict(o, boxes=F.OBSTACLES[name]) for name, o in old['obst'].items()}
    local = [L['V'][L['F']] for L in old['links']]
    path = OUT/'refined_collision_audit.json'
    results = json.loads(path.read_text()) if path.exists() else {}
    files = [Path(p) for p in args.files] if args.files else sorted(OUT.glob('*_witness.npz'))
    for source in files:
        if source.stem in results:
            continue
        data = np.load(source)
        # Other saved states already have positive coarse lower bounds.
        indices = np.flatnonzero(data['lower'] <= 0)
        rec = dict(source=source.name, saved_states=len(data['q']), inconclusive_states=len(indices),
                   checked_states=0, status='separated_at_saved_states', unresolved_states=[])
        started = time.monotonic()
        for k in indices:
            transforms, _, _ = F.fk(data['q'][k])
            hit = False
            for i, mesh in enumerate(local, 1):
                T = transforms[i]
                world = mesh @ T[:3, :3].T + T[:3, 3]
                for name, obstacle in objects.items():
                    if np.any(world.max(axis=(0, 1)) < obstacle['blo']) or np.any(world.min(axis=(0, 1)) > obstacle['bhi']):
                        continue
                    outcome = audit_pair(world, obstacle)
                    if outcome['status'] == 'collision':
                        rec.update(status='collision_witness', state=int(k), link=i, obstacle=name, witness=outcome)
                        hit = True
                        break
                    if outcome['status'] == 'unresolved':
                        rec['unresolved_states'].append(dict(state=int(k), link=i, obstacle=name, outcome=outcome))
                if hit:
                    break
            rec['checked_states'] += 1
            if hit:
                break
        if rec['status'] != 'collision_witness' and rec['unresolved_states']:
            rec['status'] = 'unresolved'
        rec['elapsed_s'] = time.monotonic()-started
        results[source.stem] = rec
        path.write_text(json.dumps(results,indent=2),encoding='utf-8')
        print(source.name, rec['status'], rec['checked_states'], '/',len(indices),
              f'{rec["elapsed_s"]:.1f}s', flush=True)


if __name__ == '__main__':
    main()
