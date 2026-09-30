"""
Numerical check of the planar summed-field coefficient bound (dock_cover.make_summed_fn): along the corner-barrier
docking trajectory (dock_cover_traj_5.npz), for random points x in every box of the craft cover whose lifted
barrier is small and random inputs u of both bodies, the true Psi~(x, u) must be >= min_k (A_k u + T_k|u| + C_k).

    python prototype_3d/summed_check2d.py
"""
import os

import numpy as np

import franka3d as F
import summed as SM
from dock_cover import J2, C, body_field, cover_2d, tight_level_2d
from dock_shapes import craft_body, station_body


def main():
    GS, GC = station_body(), craft_body()
    fS, fC = body_field(GS, .0125, .15), body_field(GC, .01, .06)
    MS, MC = F.hessian_majorant(fS), F.hessian_majorant(fC)
    lS, lC = tight_level_2d(fS, MS, GS), tight_level_2d(fC, MC, GC)
    PC, sz = cover_2d(fC, lC, .005)
    body = SM.prepare_body(fC, lC, PC, sz, 2)
    M3S = SM.field_bound(fS, sz * np.sqrt(2) / 2)
    traj = np.load(os.path.join(F.HERE, 'dock_cover_traj_5.npz'))['traj']
    rng = np.random.default_rng(1)
    Wg = np.zeros((6, 2, 2)); Wg[2], Wg[5] = -J2, J2
    worst, n = np.inf, 0
    for k in range(0, len(traj), 40):
        xi, xj = traj[k, 0], traj[k, 1]
        xi = xi + rng.normal(0, [.0, .0, .0])
        Ri, Rj = C.rot(xi[2]), C.rot(xj[2])
        cw = PC @ Rj.T + xj[:2]
        cb = (cw - xi[:2]) @ Ri
        r = body['r']
        dlo, dhi = fS.lo[:2], (fS.lo + fS.K * fS.h)[:2]
        ok = np.flatnonzero(np.all((cb - r > dlo) & (cb + r < dhi), axis=1))
        if not len(ok):
            continue
        ok = ok[fS.eval(np.c_[cb[ok], np.zeros(len(ok))]) < .03]
        if not len(ok):
            continue
        Jc = np.zeros((len(ok), 6, 2))
        Jc[:, 0], Jc[:, 1] = -Ri.T[:, 0], -Ri.T[:, 1]
        Jc[:, 3], Jc[:, 4] = Ri.T[:, 0], Ri.T[:, 1]
        Jc[:, 2] = -((cw[ok] - xi[:2]) @ J2.T) @ Ri
        Jc[:, 5] = ((cw[ok] - xj[:2]) @ J2.T) @ Ri
        Jp = J2 @ (xj[:2] - xi[:2])
        E = np.array([[-1., 0., -Jp[0], 1., 0., 0.], [0., -1., -Jp[1], 0., 1., 0.], [0., 0., -1., 0., 0., 1.]])
        Sb = np.c_[np.ones((len(ok), 2)), np.linalg.norm(cw[ok] - xj[:2], axis=1) + r]
        res = SM.rows(body, ok, Ri.T @ Rj, cb[ok], fS, lS, M3S, Jc, Wg, Sb, 1e9, np.ones(6), np.ones(3),
                      prune_rows=False, depth=0)
        if res is None:
            continue
        A, T, Cc, _, _, kept = res
        K = {'bernstein': 9, 'vertex': 4}.get(SM.ROW_MODE, 1)
        A, T, Cc = A.reshape(-1, K, 6), T.reshape(-1, K, 3), Cc.reshape(-1, K)
        for bi, b in enumerate(ok[kept]):
            u = rng.uniform(-1, 1, (30, 6)) * [1, 1, 2, 1, 1, 2]
            a = np.abs(u @ E.T)
            lower = (np.einsum('kj,nj->nk', A[bi], u) + np.einsum('kj,nj->nk', T[bi], a) + Cc[bi]).min(1)
            x = PC[b] + sz * (rng.uniform(size=(30, 2)) - .5)            # craft frame
            w = x @ Rj.T + xj[:2]
            y = (w - xi[:2]) @ Ri
            v, g = fS.eval(np.c_[y, np.zeros(30)], order=1)
            wd = u[:, 3:5] + u[:, 5:6] * ((w - xj[:2]) @ J2.T) - u[:, 0:2] - u[:, 2:3] * ((w - xi[:2]) @ J2.T)
            yd = wd @ Ri
            psi = np.einsum('pa,pa->p', g[:, :2], yd) + 5. * (fC.eval(np.c_[x, np.zeros(30)]) - lC + v - lS)
            worst = min(worst, float((psi - lower).min()))
            n += 30
    print(f'checked {n} (x, u) samples; min(Psi_true - bound) = {worst:.3e}  (must be >= 0)')


if __name__ == '__main__':
    main()
