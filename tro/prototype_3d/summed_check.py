"""
Numerical check of the summed-field coefficient bound (summed.py): at configurations of the star-tube run,
for random points x in every active box and random admissible inputs u, the true
    Psi~(x, u) = grad phi_B(y(x)) . ydot(x, u) + gamma (phi_A(x) - l_A + phi_B(y(x)) - l_B)
must be >= min_k (A_k u + T_k |u| + C_k) over the box's Bernstein coefficient rows.

    python prototype_3d/summed_check.py
"""
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'interactive'))
import server as SV                                     # noqa: E402
import franka3d as F                                    # noqa: E402
import startube_setup as ST                             # noqa: E402
import summed as S                                      # noqa: E402
import summed3d as S3                                   # noqa: E402


def main():
    path = sorted(glob.glob(os.path.join(HERE, 'interactive', 'setups', '*.json')))[-1]
    d, q0, pB, RB = ST.load(path)
    sim = SV.Sim()
    links, B = sim.links, sim.objlib['star tube']
    for L in links:
        S3.prep_link(L)
    S3.prep_field(B)
    qs = np.load(os.path.join(HERE, 'startube_setup_ours.npz'))['q']
    rng = np.random.default_rng(0)
    worst, n_checked = np.inf, 0
    for k in range(0, len(qs), 100):
        q = qs[k]
        T1, Z1, O1 = F.fk(q)
        for A in links:
            RA, pA = T1[A['frame']][:3, :3], T1[A['frame']][:3, 3]
            RbA, tb = RB.T @ RA, RB.T @ (pA - pB)
            C = A['PC'] @ RbA.T + tb
            r = A['sb']['r']
            ok = np.flatnonzero(np.all((C - r > B['dlo']) & (C + r < B['dhi']), axis=1))
            if not len(ok):
                continue
            ok = ok[B['fs'].eval(C[ok]) < .04]
            if not len(ok):
                continue
            nj = A['n_joints']
            cw = C[ok] @ RB.T + pB
            Jw = np.zeros((len(ok), 7, 3)); ww = np.zeros((7, 3))
            Jw[:, :nj] = np.cross(Z1[None, :nj], cw[:, None] - O1[None, :nj]); ww[:nj] = Z1[:nj]
            E = np.zeros((6, 7))
            E[:3, :nj] = np.cross(Z1[:nj], pA - O1[:nj]).T
            E[3:, :nj] = Z1[:nj].T
            Sb = np.c_[np.ones((len(ok), 3)), np.repeat((np.linalg.norm(cw - pA, axis=1) + r)[:, None], 3, axis=1)]
            res = S.rows(A['sb'], ok, RbA, C[ok], B['fs'], B['l'], B['M3'], Jw @ RB, S3.skew(ww @ RB), Sb, 1e9,
                         np.ones(7), np.ones(6), prune_rows=False, depth=0)
            if res is None:
                continue
            Ar, Tr, Cr, _, _, kept = res
            K = {'bernstein': 27, 'vertex': 8}.get(S.ROW_MODE, 1)
            Ar, Tr, Cr = Ar.reshape(-1, K, 7), Tr.reshape(-1, K, 6), Cr.reshape(-1, K)
            for bi, b in enumerate(ok[kept]):
                u = rng.uniform(-1, 1, (20, 7)); u[:, nj:] = 0
                a = np.abs(u @ E.T)
                lower = (np.einsum('kj,nj->nk', Ar[bi], u) + np.einsum('kj,nj->nk', Tr[bi], a) + Cr[bi]).min(1)
                x = A['PC'][b] + A['side'] * (rng.uniform(size=(20, 3)) - .5)          # link frame
                y = x @ RbA.T + tb
                w = x @ RA.T + pA
                v, g = B['fs'].eval(y, order=1)
                Jx = np.cross(Z1[None, :nj], w[:, None] - O1[None, :nj])            # (20, nj, 3) world
                yd = np.einsum('pja,pj->pa', Jx, u[:, :nj]) @ RB                      # tube frame
                psi = np.einsum('pa,pa->p', g, yd) + SV.GAMMA * (A['fs_'].eval(x) - A['l'] + v - B['l'])
                worst = min(worst, float((psi - lower).min()))
                n_checked += 20
    print(f'checked {n_checked} (x, u) samples; min(Psi_true - bound) = {worst:.3e}  (must be >= 0)')


if __name__ == '__main__':
    main()
