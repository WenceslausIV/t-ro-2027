"""
Timing study: the star-tube example (startube_setup.py) with the SUM-OF-SDF barrier family and Bernstein
coefficient constraints on Psi = hdot + gamma h (safe set continuous; boxes only partition the parameter domain).

For every x in the shell of link A (covered by A's boxes, in A's frame),
    h(x) = phi_A(x) - l_A + phi_B(y(x)) - l_B,   y(x) = link point x in the tube frame,
and per active box the 64 tricubic Bernstein coefficients of a lower model of Psi(x, qd) must be >= 0:
phi_A exactly (tricubic on the box), phi_B by its second-order Taylor model at the box center with a certified
third-derivative remainder (value: M3 r^3 / 6, gradient term: M3 r^2 / 2 times the maximal point speed).
Each coefficient is affine in the joint velocities. Boxes are skipped when a certified lower bound of h on them
is >= eta (the activation argument of the paper); rows that hold for every admissible input are dropped.

    python prototype_3d/star_bern_timing.py            -> star_bern_timing.json
"""
import glob
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'interactive'))
import server as S                                      # noqa: E402
import franka3d as F                                    # noqa: E402
import startube_setup as ST                             # noqa: E402
from proto3d import Spline3                             # noqa: E402
from sdf_cbf_utils import solve_ldp_qp                  # noqa: E402

EPS = 1e-3                                              # weight of the auxiliary speed bounds t in the QP
XI = np.array([0., 1 / 3, 2 / 3, 1.])
TT = XI - .5                                            # Bernstein coefficients of (xi - 1/2) in the cubic basis
SQ = np.array([.25, -1 / 12, -1 / 12, .25])             # ... of (xi - 1/2)^2
K3 = np.stack(np.meshgrid(np.arange(4), np.arange(4), np.arange(4), indexing='ij'), -1).reshape(-1, 3)
TK = TT[K3]                                             # (64, 3)
QT = TK[:, :, None] * TK[:, None, :]                    # (64, 3, 3): coefficients of t_a t_b (a != b)
for a in range(3):
    QT[:, a, a] = SQ[K3[:, a]]                          # ... and of t_a^2
BM = np.array([[(1 - x) ** 3, 3 * x * (1 - x) ** 2, 3 * x ** 2 * (1 - x), x ** 3] for x in XI])
BMI = np.linalg.inv(BM)                                 # values at XI -> cubic Bernstein coefficients


def third_majorant(f):
    """C^2 majorant of a certified cellwise bound on the Frobenius norm of the third-derivative tensor."""
    C = f.cell_coeffs()
    d = lambda X, a: (X.shape[3 + a] - 1) * np.diff(X, axis=3 + a) / f.h
    s = 0.
    for a in range(3):
        for b in range(a, 3):
            for c in range(b, 3):
                m = np.abs(d(d(d(C, a), b), c)).max(axis=(3, 4, 5))
                mult = 6 / (1 + (a == b) + (b == c) + (a == c) + 2 * (a == b == c))   # ordered permutations
                s = s + mult * m ** 2
    W = np.pad(np.sqrt(s), 4, mode='edge')
    for a in range(3):
        W = np.lib.stride_tricks.sliding_window_view(W, 6, axis=a).max(axis=-1)
    g = Spline3.__new__(Spline3)
    g.lo, g.h, g.K, g.n, g.W = f.lo, f.h, f.K, f.n, W
    return g


def prepare(links):
    for A in links:
        s = A['side']
        lo = A['PC'] - s / 2
        P = (lo[:, None, :] + s * XI[K3][None]).reshape(-1, 3)
        V = A['fs'].eval(P).reshape(-1, 4, 4, 4)
        A['BC'] = np.einsum('ip,jq,kr,npqr->nijk', BMI, BMI, BMI, V).reshape(-1, 64) - A['l']
        A['bmin'] = A['BC'].min(axis=1)
        A['gmin'] = np.array([A['bmin'][g].min() for g in A['groups']])
        A['goff'] = s * TK                              # Greville offsets from the box center (link frame)


def bern_rows(A, TA, Z1, O1, B, RB, pB, eta, qd_max):
    RA, pA = TA[:3, :3], TA[:3, 3]
    if np.linalg.norm(pA - pB) > A['rad'] + B['reach'] + eta:
        return None
    RbA, tb = RB.T @ RA, RB.T @ (pA - pB)
    s, r = A['side'], A['r']
    cc = A['cc'] @ RbA.T + tb
    keep = ~np.all((cc - A['crad'][:, None] > B['dlo']) & (cc + A['crad'][:, None] < B['dhi']), axis=1)
    ins = np.flatnonzero(~keep)
    if len(ins):
        ci = B['fs'].cell_of(cc[ins])
        low = (B['fs'].eval(cc[ins]) - B['Gn'][ci[:, 0], ci[:, 1], ci[:, 2]] * A['crad'][ins]
               - .5 * B['Mn'][ci[:, 0], ci[:, 1], ci[:, 2]] * r ** 2 - B['l'] + A['gmin'][ins])
        keep[ins[low < eta]] = True
    if not keep.any():
        return None
    sub = np.concatenate([A['groups'][j] for j in np.flatnonzero(keep)])
    C = A['PC'][sub] @ RbA.T + tb
    ok = np.all((C - r > B['dlo']) & (C + r < B['dhi']), axis=1)
    idx = sub[ok]
    C = C[ok]
    if not len(idx):
        return None
    val, g = B['fs'].eval(C, order=1)
    pre = val - np.linalg.norm(g, axis=1) * r - .5 * B['Ms'].eval(C) * r ** 2 - B['l'] + A['bmin'][idx] < eta
    idx, c = idx[pre], C[pre]
    if not len(idx):
        return None
    v, g, H = B['fs'].eval(c, order=2)
    M3 = B['M3'].eval(c)
    Hl = np.einsum('ai,nab,bj->nij', RbA, H, RbA)                   # Hessian in the box (link) axes
    d = A['goff'] @ RbA.T                                           # (64, 3) Greville offsets, tube frame
    hk = (A['BC'][idx] + (v - B['l'] - M3 * r ** 3 / 6)[:, None] + g @ d.T
          + .5 * s ** 2 * np.einsum('kab,nab->nk', QT, Hl))         # Bernstein coefficients of a lower model of h
    act = hk.min(axis=1) < eta
    if not act.any():
        return None
    idx, c, v, g, H, M3, hk = idx[act], c[act], v[act], g[act], H[act], M3[act], hk[act]
    nj = A['n_joints']
    cw = c @ RB.T + pB                                              # box centers, world
    om = Z1[:nj] @ RB                                               # (nj, 3) joint axes, tube frame
    Jc = np.einsum('ab,njb->nja', RB.T, np.cross(Z1[None, :nj], cw[:, None] - O1[None, :nj]))  # (n, nj, 3)
    # gradient-model velocity term: (g + H d)^T (Jc + om x d)
    Hd = np.einsum('nab,kb->nka', H, d)                             # (n, 64, 3)
    A1 = np.einsum('nka,nja->nkj', g[:, None] + Hd, Jc)
    A1 += np.einsum('ja,nka->nkj', om, np.cross(d[None], g[:, None]))          # g.(om x d) = om.(d x g)
    # quadratic part d^T H [om]x d, expressed through its Bernstein coefficients in the box axes
    Wx = np.zeros((nj, 3, 3))
    Wx[:, 0, 1], Wx[:, 0, 2], Wx[:, 1, 2] = -om[:, 2], om[:, 1], -om[:, 0]
    Wx -= Wx.transpose(0, 2, 1)
    Mq = np.einsum('nab,jbc->njac', H, Wx)
    Mq = np.einsum('ai,njab,bk->njik', RbA, Mq, RbA)
    A1 += s ** 2 * np.einsum('kab,njab->nkj', QT, Mq)
    # gradient remainder: |R2^T ydot(x)| <= (M3 r^2 / 2) sum_j Lj |qd_j| <= (M3 r^2 / 2) sum_j Lj t_j, t_j >= |qd_j|,
    # Lj >= max over the box of the point-velocity column norm; vanishes at rest
    Lj = np.linalg.norm(Jc, axis=2) + np.linalg.norm(om, axis=1)[None] * r          # (n, nj)
    Tt = np.repeat(-.5 * (M3 * r ** 2)[:, None] * Lj, 64, axis=0)                     # (n*64, nj)
    C0 = (S.GAMMA * hk).ravel()
    A1 = A1.reshape(-1, nj)
    need = C0 < (np.abs(A1) - Tt) @ qd_max[:nj]                     # rows that some admissible input can violate
    rows = np.zeros((need.sum(), 14))
    rows[:, :nj], rows[:, 7:7 + nj] = A1[need], Tt[need]
    return rows, -C0[need], hk.min(), len(idx)


def main():
    path = sorted(glob.glob(os.path.join(HERE, 'interactive', 'setups', '*.json')))[-1]
    d, q0, pB, RB = ST.load(path)
    sim = S.Sim()
    links, B = sim.links, sim.objlib['star tube']
    t0 = time.perf_counter()
    B['M3'] = third_majorant(B['fs'])
    prepare(links)
    print(f'prepare {time.perf_counter() - t0:.1f} s; boxes per link', [len(A['PC']) for A in links], flush=True)
    Th, _, _ = F.fk(S.Q_HOME)
    goal_p, goal_R = Th[8][:3, :3] @ S.TCP + Th[8][:3, 3], Th[8][:3, :3]
    eta = F.ACT
    qd_max = np.broadcast_to(S.QD_MAX, 7).astype(float)
    out = dict(setup=os.path.basename(path), eta=eta, gamma=S.GAMMA, dt=ST.DT)
    for method in ('cover', 'bern'):
        q, u, z_prev = q0.copy(), None, None
        L = dict(t_rows=[], t_qp=[], rows=[], boxes=[], gap=[], hmin=[], slack=0)
        for k in range(ST.STEPS):
            T1, Z1, O1 = F.fk(q)
            p_ee = T1[8][:3, :3] @ S.TCP + T1[8][:3, 3]
            J = np.vstack([np.cross(Z1, p_ee - O1).T, Z1.T])
            e = np.r_[np.clip(2. * (goal_p - p_ee), -.5, .5), np.clip(1.5 * S._log(goal_R @ T1[8][:3, :3].T), -1., 1.)]
            Jp = J.T @ np.linalg.inv(J @ J.T + 1e-3 * np.eye(6))
            u_nom = np.clip(Jp @ e + (np.eye(7) - Jp @ J) @ (.5 * (S.Q_HOME - q)), -S.QD_MAX, S.QD_MAX)
            t0 = time.perf_counter()
            rows, rhs, nb, hm = [], [], 0, np.inf
            if method == 'cover':
                Z1xO1 = np.cross(Z1, O1)
                for A in links:
                    res = S.body_rows(A, T1[A['frame']], Z1, Z1xO1, B, RB, pB)
                    if res is not None:
                        rows.append(res[0]); rhs.append(-S.GAMMA * res[1]); hm = min(hm, res[1].min())
            else:
                for A in links:
                    res = bern_rows(A, T1[A['frame']], Z1, O1, B, RB, pB, eta, qd_max)
                    if res is not None:
                        rows.append(res[0]); rhs.append(res[1]); hm = min(hm, res[2]); nb += res[3]
            t1 = time.perf_counter()
            lo = np.clip(-2. * (q - F.Q_MIN), -S.QD_MAX, 0.)
            hi = np.clip(2. * (F.Q_MAX - q), 0., S.QD_MAX)
            n_bar = int(sum(len(r) for r in rows))
            if method == 'cover':
                G = np.vstack(rows + [np.eye(7), -np.eye(7)])
                hv = np.r_[np.concatenate(rhs) if rhs else np.zeros(0), lo, -hi]
                u, sl = solve_ldp_qp(u_nom, G, hv, n_bar, u)
            else:
                # z = (qd, sqrt(EPS) t): min |qd - qd_nom|^2 + EPS |t|^2, t >= |qd| shared by all boxes
                se = np.sqrt(EPS)
                I, Z = np.eye(7), np.zeros((7, 7))
                Gb = np.vstack(rows) if rows else np.zeros((0, 14))
                Gb = np.c_[Gb[:, :7], Gb[:, 7:] / se]
                G = np.vstack([Gb, np.c_[-I, I / se], np.c_[I, I / se], np.c_[I, Z], np.c_[-I, Z]])
                hv = np.r_[np.concatenate(rhs) if rhs else np.zeros(0), np.zeros(14), lo, -hi]
                z, sl = solve_ldp_qp(np.r_[u_nom, np.zeros(7)], G, hv, n_bar, z_prev)
                z_prev, u = z, z[:7]
            t2 = time.perf_counter()
            L['t_rows'].append(t1 - t0); L['t_qp'].append(t2 - t1); L['rows'].append(n_bar); L['boxes'].append(nb)
            L['hmin'].append(hm); L['slack'] += sl > 0
            L['gap'].append(ST.gap(q, links, B, RB, pB))
            q = q + ST.DT * u
            if k % 250 == 0:
                print(f'  [{method}] step {k}: rows {n_bar} boxes {nb} build {1e3 * (t1 - t0):.1f} ms qp {1e3 * (t2 - t1):.1f} ms '
                      f'gap {1e3 * L["gap"][-1]:.1f} mm', flush=True)
        tr, tq = 1e3 * np.array(L['t_rows']), 1e3 * np.array(L['t_qp'])
        tot = tr + tq
        act = np.array(L['rows']) > 0
        r = dict(min_gap_mm=1e3 * min(L['gap']), slack_steps=int(L['slack']), rows_max=int(max(L['rows'])),
                 rows_med_active=float(np.median(np.array(L['rows'])[act])) if act.any() else 0,
                 boxes_max=int(max(L['boxes'])),
                 total_med_ms=float(np.median(tot)), total_p95_ms=float(np.percentile(tot, 95)), total_max_ms=float(tot.max()),
                 active_med_ms=float(np.median(tot[act])) if act.any() else 0,
                 active_p95_ms=float(np.percentile(tot[act], 95)) if act.any() else 0,
                 build_med_ms=float(np.median(tr[act])) if act.any() else 0,
                 qp_med_ms=float(np.median(tq[act])) if act.any() else 0,
                 qp_p95_ms=float(np.percentile(tq[act], 95)) if act.any() else 0,
                 steps_active=int(act.sum()), final_q=q.tolist())
        out[method] = r
        print(method, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items() if k != 'final_q'}, flush=True)
        np.savez_compressed(os.path.join(HERE, f'star_bern_{method}.npz'), rows=np.array(L['rows']), t_rows=tr, t_qp=tq,
                            gap=np.array(L['gap']), boxes=np.array(L['boxes']))
    json.dump(out, open(os.path.join(HERE, 'star_bern_timing.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
