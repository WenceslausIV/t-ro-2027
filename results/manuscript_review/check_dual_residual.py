"""Numerical falsification check of two-moving-arm derivatives and finite bounds.

Uses existing cached fields and interpolated historical start/goal pairs. This is
not a trajectory safety audit or a proof by sampling. Original artifacts are read only.
"""
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
import dual3d as D
import summed as S
import summed3d as S3


def main():
    links, _, _ = F.build(cover_side='cached')
    for L in links:
        L['fs'] = F.spline_from(L['f'])
        L['Ms'] = F.hessian_majorant(L['fs'])
        L['dlo'], L['dhi'] = L['fs'].lo, L['fs'].lo+L['fs'].K*L['fs'].h
    reach = max(L['crad'].max() for L in links)
    for L in links:
        L['Gn'], _ = F.neighborhood_bounds(L['fs'], reach)
    S3.prep_all(links)
    trials = json.loads((ROOT/'prototype_3d/dual_trials.json').read_text())
    rng = np.random.default_rng(1943)
    modes = [(r, w) for r in ('vertex', 'bernstein') for w in ('one', 'norm', 'proj')]
    minima = {f'{r}-{w}': np.inf for r, w in modes}
    velocity_error, derivative_error, n = 0., 0., 0
    pairs = []
    original_rows = S.rows
    captured = []

    def capture(*args, **kwargs):
        a = list(args)
        take = np.arange(min(8, len(a[1])))
        for k in (1, 3, 7, 9):
            a[k] = a[k][take]
        a[10] = 1e12
        captured.append(a)
        return original_rows(*a, prune_rows=False, depth=0, row_mode='vertex', mult_mode='one')

    for trial in (0, 2, 7):
        q0, qg = np.asarray(trials[trial])
        for fraction in (.35, .65):
            q = (1-fraction)*q0+fraction*qg
            (T1, Z1, O1), (T2, Z2, O2) = D.fk2(q)
            done = 0
            for A in links[::-1]:
                if done >= 3:
                    break
                RA, pA = T1[A['frame']][:3,:3], T1[A['frame']][:3,3]
                for B in links[::-1]:
                    RB, pB = T2[B['frame']][:3,:3], T2[B['frame']][:3,3]
                    captured.clear()
                    S.rows = capture
                    try:
                        result = S3.pair_rows(A, RA, pA, B, RB, pB,
                            [(Z1,O1,A['n_joints'],1.,0),(Z2,O2,B['n_joints'],-1.,7)],
                            F.ACT, np.ones(14))
                    finally:
                        S.rows = original_rows
                    if result is None or not captured:
                        continue
                    a, E = captured[0], result[-1]
                    body, ids, R, yc, fB, level = a[:6]
                    count = len(ids)
                    # One independently perturbed material point/input per retained box.
                    for _ in range(6):
                        x = body['PC'][ids] + body['side']*(rng.random((count,3))-.5)
                        y = x@R.T + RB.T@(pA-pB)
                        u = rng.uniform(-1,1,14)
                        yd = np.einsum('nmd,m->nd',a[7],u) + (y-yc)@np.einsum('m,mij->ij',u,a[8]).T
                        dt = 1e-6
                        evaluated = []
                        for sign in (-1,1):
                            (Ta,_,_),(Tb,_,_) = D.fk2(q+sign*dt*u)
                            ta, tb = Ta[A['frame']], Tb[B['frame']]
                            evaluated.append((x@ta[:3,:3].T+ta[:3,3]-tb[:3,3])@tb[:3,:3])
                        fd = (evaluated[1]-evaluated[0])/(2*dt)
                        velocity_error = max(velocity_error,float(np.abs(yd-fd).max()))
                        v, g = fB.eval(y,order=1)
                        hdot = np.einsum('nd,nd->n',g,yd)
                        fd_h = (fB.eval(evaluated[1])-fB.eval(evaluated[0]))/(2*dt)
                        derivative_error = max(derivative_error,float(np.abs(hdot-fd_h).max()))
                        ga = A['fs_'].eval(x)-A['l']
                        _, center_g = fB.eval(yc,order=1)
                        for row, mult in modes:
                            ar,tr,cr,_,_,kept = original_rows(*a,prune_rows=False,depth=0,row_mode=row,mult_mode=mult)
                            assert len(kept)==count
                            lower = (ar@u+tr@np.abs(E@u)+cr).reshape(count,-1).min(axis=1)
                            weight = S.multiplier_weights(body['gA'][ids],center_g@R,mult)
                            residual = hdot+F.GAMMA*(weight*ga+v-level)
                            minima[f'{row}-{mult}'] = min(minima[f'{row}-{mult}'],float((residual-lower).min()))
                        n += count
                    pairs.append([trial,fraction,A['frame'],B['frame'],count])
                    done += 1
                    if done >= 3:
                        break
    assert n > 100 and len(pairs) >= 6, (n,pairs)
    assert velocity_error < 1e-7 and derivative_error < 1e-6
    assert min(minima.values()) >= -1e-9, minima
    out = dict(scope=__doc__,point_input_samples=n,pairs=pairs,
               max_velocity_finite_difference_error=velocity_error,
               max_field_derivative_finite_difference_error=derivative_error,
               minimum_true_residual_minus_bound=minima)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k not in ('scope','pairs')},indent=2))


if __name__ == '__main__':
    main()
