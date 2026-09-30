"""
Summed-field barrier rows (summed.py) for the spatial experiments: a link of an arm against a static object
(Franka obstacles, star tube) or against a link of a second arm (dual arm). Pruning follows the paper:
clusters and boxes are skipped when a certified lower bound of the lifted barrier h~ on them is >= eta.
"""
import numpy as np

import franka3d as F
import summed as S


def _get(B, *keys):
    for k in keys:
        if k in B:
            return B[k]
    raise KeyError(keys)


def prep_link(L):
    """Offline data of a link (A role): Taylor data on its boxes and per-cluster lower bounds of phi_A - l_A."""
    if 'sb' in L:
        return
    f = L['fs'] if 'fs' in L else F.spline_from(L['f'])
    L['fs_'] = f
    L['sb'] = S.prepare_body(f, L['l'], L['PC'], L['side'], 3)
    L['gminA'] = np.array([L['sb']['lowA'][g].min() for g in L['groups']])


def prep_field(B, reach=None):
    """Build a radius-specific majorant; defer until pair_rows when reach is not yet known."""
    if reach is None:
        return
    if 'M3' not in B or B.get('_M3_reach', -1.) < reach:
        B['M3'] = S.field_bound(_get(B, 'fs', 'f'), reach)
        B['_M3_reach'] = reach
    f = _get(B, 'fs', 'f')
    if reach > f.h and B.get('_M2_reach', -1.) < reach:
        B['M2_radius'] = S.LinearMajorant(f, f.hessian_bound(cellwise=True), reach)
        B['_M2_reach'] = reach


def skew(w):
    W = np.zeros(w.shape[:-1] + (3, 3))
    W[..., 0, 1], W[..., 0, 2], W[..., 1, 2] = -w[..., 2], w[..., 1], -w[..., 0]
    return W - np.swapaxes(W, -1, -2)


def pair_rows(A, RA, pA, B, RB, pB, joints, eta, umax, gamma=F.GAMMA, sd=None):
    """Rows of link A (world pose RA, pA) against the field of B (world pose RB, pB).
    joints: list of (Z (k,3), O (k,3), n_moving, sign, column offset) describing which inputs move A (+) and B (-).
    Returns (A_rows, T_rows, C, n_boxes, h_low) or None."""
    sb, r = A['sb'], A['sb']['r']
    prep_field(B, r)
    fB, lB, MsB, M3B = _get(B, 'fs', 'f'), B['l'], _get(B, 'Ms', 'M'), B['M3']
    if r > fB.h:
        MsB = B['M2_radius']
    RbA, tb = RB.T @ RA, RB.T @ (pA - pB)
    cc = A['cc'] @ RbA.T + tb
    inside = np.all((cc - A['crad'][:, None] > B['dlo']) & (cc + A['crad'][:, None] < B['dhi']), axis=1)
    keep = ~inside                                  # clusters leaving the domain are checked box by box
    ins = np.flatnonzero(inside)
    if len(ins):
        ci = fB.cell_of(cc[ins])
        low = fB.eval(cc[ins]) - B['Gn'][ci[:, 0], ci[:, 1], ci[:, 2]] * A['crad'][ins] - lB + A['gminA'][ins]
        keep[ins[low < eta]] = True
    if not keep.any():
        return None
    sub = np.concatenate([A['groups'][j] for j in np.flatnonzero(keep)])
    C = A['PC'][sub] @ RbA.T + tb
    ok = np.all((C - r > B['dlo']) & (C + r < B['dhi']), axis=1)
    sub, C = sub[ok], C[ok]
    if not len(sub):
        return None
    val, g = fB.eval(C, order=1)
    pre = val - np.linalg.norm(g, axis=1) * r - .5 * MsB.eval(C) * r ** 2 - lB + sb['lowA'][sub] < eta
    sub, C = sub[pre], C[pre]
    if not len(sub):
        return None
    cw = C @ RB.T + pB
    m = max(j[4] + len(j[0]) for j in joints)
    Jw = np.zeros((len(sub), m, 3))
    ww = np.zeros((m, 3))
    for Z, O, nmov, sign, off in joints:
        Jw[:, off:off + nmov] = sign * np.cross(Z[None, :nmov], cw[:, None] - O[None, :nmov])
        ww[off:off + nmov] = sign * Z[:nmov]
    Jc, W = Jw @ RB, skew(ww @ RB)
    # relative twist of A with respect to B at A's origin (world): the speed of every box point is at most
    # |V| + |Omega| (|c - p_A| + r) <= (sum_a a_a) + (|c - p_A| + r) (sum_a a_{3+a}),  a >= |E u|
    E = np.zeros((6, m))
    for Z, O, nmov, sign, off in joints:
        E[:3, off:off + nmov] = sign * np.cross(Z[:nmov], pA - O[:nmov]).T
        E[3:, off:off + nmov] = sign * Z[:nmov].T
    rho = np.linalg.norm(cw - pA, axis=1) + r
    Sb = np.c_[np.ones((len(sub), 3)), np.repeat(rho[:, None], 3, axis=1)]
    if sd is not None:                              # sampled data: auxiliary (a, s) with s >= |u|, and caps
        sd = dict(sd, M2=B['M2_sd'], m=m)
        E = np.r_[E, np.eye(m)]
    res = S.rows(sb, sub, RbA, C, fB, lB, M3B, Jc, W, Sb, eta, umax, np.abs(E) @ umax, gamma, sd=sd)
    if res is None:
        return None
    if sd is not None:                              # 1.a_V <= V, 1.a_Omega <= Om, 1.s <= nu
        k = len(E)
        Tcap = np.zeros((3, k))
        Tcap[0, :3], Tcap[1, 3:6], Tcap[2, 6:] = -1., -1., -1.
        res = list(res)
        res[0] = np.r_[res[0], np.zeros((3, m))]
        res[1] = np.r_[res[1], Tcap]
        res[2] = np.r_[res[2], [sd['V'], sd['Om'], sd['nu']]]
        if len(res) == 7:
            import scipy.sparse as sp
            res[6] = sp.vstack([res[6], sp.csr_matrix((3, res[6].shape[1]))], format='csr')
        res = tuple(res)
    Ar, Tr, Cr, nbox, hlow = res[:5]
    if len(res) == 7:                               # 'free' multipliers: one column per active box
        return Ar, Tr, Cr, nbox, hlow, E, res[6]
    return Ar, Tr, Cr, nbox, hlow, E


def stack(results, m):
    """Stack the rows of several pairs; each pair has its own block of 6 auxiliary variables."""
    res = [r for r in results if r is not None]
    if not res:
        return np.zeros((0, m)), np.zeros((0, 0)), np.zeros(0), 0, np.inf, np.zeros((0, m))
    k = sum(len(r[5]) for r in res)
    if len(res[0]) == 7:                            # 'free' multipliers: block-diagonal T kept sparse
        import scipy.sparse as sp
        T = sp.block_diag([sp.csr_matrix(r[1]) for r in res], format='csr')
    else:
        T = np.zeros((sum(len(r[2]) for r in res), k))
        i = j = 0
        for r in res:
            T[i:i + len(r[2]), j:j + len(r[5])] = r[1]
            i += len(r[2]); j += len(r[5])
    out = (np.vstack([r[0] for r in res]), T, np.concatenate([r[2] for r in res]),
           sum(r[3] for r in res), min(r[4] for r in res), np.vstack([r[5] for r in res]))
    if len(res[0]) == 7:                            # block-diagonal multiplier columns of all pairs (sparse)
        import scipy.sparse as sp
        return out + (sp.block_diag([sp.csr_matrix(r[6]) for r in res], format='csr'),)
    return out


# ---------------------------------------------------------------------------------------------
# experiment-specific row builders
# ---------------------------------------------------------------------------------------------
UMAX7 = np.full(F.DOF, F.QD_MAX)
UMAX14 = np.full(2 * F.DOF, F.QD_MAX)
I3 = np.eye(3)


def franka_rows(q, links, obst, eta=F.ACT):
    """All (link, obstacle) pairs of the Franka experiment; obstacles are static, in the world frame."""
    T, Z, O = F.fk(q)
    res = []
    for L in links:
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        for o in obst.values():
            if np.any(p + L['rad'] < o['dlo']) or np.any(p - L['rad'] > o['dhi']):
                continue
            res.append(pair_rows(L, R, p, o, I3, np.zeros(3), [(Z, O, L['n_joints'], 1., 0)], eta, UMAX7))
    return stack(res, F.DOF)


def dual_rows(q, links, eta=F.ACT):
    """All (link of arm 1, link of arm 2) pairs; both move, 14 inputs."""
    import dual3d as D
    (T1, Z1, O1), (T2, Z2, O2) = D.fk2(q)
    res = []
    for A in links:
        RA, pA = T1[A['frame']][:3, :3], T1[A['frame']][:3, 3]
        for B in links:
            RB, pB = T2[B['frame']][:3, :3], T2[B['frame']][:3, 3]
            if np.linalg.norm(pA - pB) > A['rad'] + B['rad'] + .06 + eta:
                continue
            res.append(pair_rows(A, RA, pA, B, RB, pB, [(Z1, O1, A['n_joints'], 1., 0),
                                                        (Z2, O2, B['n_joints'], -1., F.DOF)], eta, UMAX14))
    return stack(res, 2 * F.DOF)


def prep_all(links, fields=()):
    reach = max(L['side'] * np.sqrt(3) / 2 for L in links)
    for L in links:
        prep_link(L)
        if 'fs' in L:
            prep_field(L, reach)
    for B in fields:
        prep_field(B, reach)
