"""
Consistency check of the online refinement (summed._children): the half-size boxes carry the same data as if
they had been built directly (centers in B's frame, center velocities, Taylor data of phi_A), and every surface
point of a parent box lies in a kept child (checked on points of S_A found by bisection along random segments).

    python prototype_3d/summed_check_children.py
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
    S3.prep_all(links, [B])
    q = np.load(os.path.join(HERE, 'startube_setup_ours.npz'))['q'][250]
    T1, Z1, O1 = F.fk(q)
    rng = np.random.default_rng(0)
    A = links[-1]
    RA, pA = T1[A['frame']][:3, :3], T1[A['frame']][:3, 3]
    RbA, tb = RB.T @ RA, RB.T @ (pA - pB)
    sb, nj = A['sb'], A['n_joints']
    idx = rng.choice(len(A['PC']), 200, replace=False)
    c = A['PC'][idx]
    cw = c @ RA.T + pA
    Jw = np.zeros((len(idx), 7, 3)); ww = np.zeros((7, 3))
    Jw[:, :nj] = np.cross(Z1[None, :nj], cw[:, None] - O1[None, :nj]); ww[:nj] = Z1[:nj]
    W = S3.skew(ww @ RB)
    D = dict(c=c, yc=c @ RbA.T + tb, vA=sb['vA'][idx], gA=sb['gA'][idx], HA=sb['HA'][idx], m3A=sb['m3A'][idx],
             Jc=Jw @ RB, Sb=np.ones((len(idx), 6)))
    Ch = S._children(D, sb, sb['side'], sb['r'], RbA, W, 3)
    yc_direct = Ch['c'] @ RbA.T + tb
    cwc = Ch['c'] @ RA.T + pA
    Jd = np.zeros((len(cwc), 7, 3)); Jd[:, :nj] = np.cross(Z1[None, :nj], cwc[:, None] - O1[None, :nj])
    print('children kept', len(Ch['c']), 'of', 8 * len(idx))
    print('max |yc - direct|', float(np.abs(Ch['yc'] - yc_direct).max()),
          ' max |Jc - direct|', float(np.abs(Ch['Jc'] - Jd @ RB).max()))
    # surface points of the parent boxes lie in kept children
    f, l, s = sb['f'], sb['level'], sb['side']
    miss, tot = 0, 0
    for k in range(len(idx)):
        a = c[k] + s * (rng.uniform(size=(64, 3)) - .5)
        b = c[k] + s * (rng.uniform(size=(64, 3)) - .5)
        fa, fb = f.eval(a) - l, f.eval(b) - l
        m = fa * fb < 0
        a, b, fa = a[m], b[m], fa[m]
        for _ in range(40):                                # bisection to a point of S_A
            mid = (a + b) / 2
            fm = f.eval(mid) - l
            left = fa * fm <= 0
            b = np.where(left[:, None], mid, b); a = np.where(left[:, None], a, mid); fa = np.where(left, fa, fm)
        p = (a + b) / 2
        mine = np.all(np.abs(Ch['c'][:, None] - p[None]) <= s / 4 + 1e-12, axis=2).any(axis=0)
        miss += int((~mine).sum()); tot += len(p)
    print(f'surface points in parent boxes: {tot}, not in a kept child: {miss}  (must be 0)')


if __name__ == '__main__':
    main()
