"""Sampled-data enforcement of the surface-cover certificate for the Franka experiment.

The input is held for one step dt. For a surface point x of link L, e(tau) = phi_O(y(x, t_k + tau)) - l_O
satisfies |e''| <= |Hess phi_O| |ydot|^2 + |grad phi_O| |yddot|. With lever-arm bound D_L (every box point and
joint origin of the chain up to L lies within D_L of every earlier joint origin), constant joint velocities give
    |ydot| <= |V| + |Omega| rbar + dt |yddot|,   |yddot| <= 1.5 D_L |u|_1^2,
and caps 1.a_V <= V, 1.a_Omega <= Om, |u|_1 <= nu (enforced in the QP, adapted from the previous input) make the
bound affine in the QP variables (summed.sampled_tightening). Rows then certify e' + gamma e >= dt/2 sup|e''|,
and with gamma dt <= 1, e >= 0 on the whole step (Theorem sampled of the paper). Skipped boxes need h >= G travel:
the activation threshold becomes max(eta, G travel), and the domain collar bound must exceed G travel.
Both travel bounds need the caps, so they are enforced for every link at every step.
"""
import numpy as np

import franka3d as F
import summed as S
import summed3d as S3

GROW, FLOOR_V, FLOOR_OM, FLOOR_NU = 1.5, .05, .1, .3     # cap adaptation per step (m/s, rad/s, rad/s)
SPEED_MAX = 3.                                          # cap on the certified point-speed bound [m/s]
CONST = {}


def _collar_bound(f, level, width):
    """Certified lower bound of phi - level on all cells within `width` of the domain boundary."""
    cmin = f.cell_coeffs().min(axis=(3, 4, 5))
    w = int(np.ceil(width / f.h))
    mask = np.zeros(cmin.shape, bool)
    for a in range(3):
        idx = [slice(None)] * 3
        idx[a] = slice(0, w); mask[tuple(idx)] = True
        idx[a] = slice(cmin.shape[a] - w, None); mask[tuple(idx)] = True
    return float(cmin[mask].min() - level)


def prepare(links, obst):
    T, Z, O = F.fk(np.zeros(F.DOF))
    r_max = max(L['r'] for L in links)
    for L in links:
        n = L['n_joints']
        chain = sum(np.linalg.norm(O[i] - O[i - 1]) for i in range(1, n))
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        rho = float(np.max(np.linalg.norm(L['PC'] @ R.T + p - O[n - 1], axis=1))) + L['r']
        L['sd_D'] = chain + rho
        L['sd_rbar'] = float(np.max(np.linalg.norm(L['PC'], axis=1))) + L['r']
    travel_max = F.DT * SPEED_MAX
    for o in obst.values():
        f = o['f']
        o['M2_sd'] = S.LinearMajorant(f, f.hessian_bound(cellwise=True), r_max + travel_max)
        o['sd_G'] = float(o['Gn'].max())
        o['sd_collar'] = _collar_bound(f, o['l'], 2 * r_max)
    CONST.update(G=max(o['sd_G'] for o in obst.values()), collar=min(o['sd_collar'] for o in obst.values()),
                 travel_max=travel_max, D={L['frame']: L['sd_D'] for L in links})
    # skipped points near the domain boundary must not reach the level within one step
    assert CONST['G'] * travel_max < CONST['collar'], CONST


class State:
    def __init__(self):
        self.u = np.zeros(F.DOF)
        self.caps = None

    def update(self, u, E, z):
        self.u = np.array(u, float)

    def last_caps(self):
        return self.caps


def _link_twist_map(q, L, fk=None):
    T, Z, O = F.fk(q) if fk is None else fk
    n, p = L['n_joints'], T[L['frame']][:3, 3]
    E = np.zeros((6, F.DOF))
    E[:3, :n] = np.cross(Z[:n], p - O[:n]).T
    E[3:, :n] = Z[:n].T
    return E


def caps_for(q, L, u_prev, E=None):
    tw = (_link_twist_map(q, L) if E is None else E) @ u_prev
    V = GROW * np.abs(tw[:3]).sum() + FLOOR_V
    Om = GROW * np.abs(tw[3:]).sum() + FLOOR_OM
    nu = GROW * np.abs(u_prev).sum() + FLOOR_NU
    c2 = 1.5 * L['sd_D']
    while V + L['sd_rbar'] * Om + F.DT * c2 * nu * nu > SPEED_MAX:
        V, Om, nu = .9 * V, .9 * Om, .9 * nu
    speed = V + L['sd_rbar'] * Om + F.DT * c2 * nu * nu
    return dict(dt=F.DT, V=V, Om=Om, nu=nu, c2=c2, travel=F.DT * speed)


def _cap_block(EL, sd):
    """Caps of a link without active boxes: 1.a_V <= V, 1.a_Omega <= Om, 1.s <= nu with a >= |E u|, s >= |u|."""
    E = np.r_[EL, np.eye(F.DOF)]
    Tc = np.zeros((3, len(E)))
    Tc[0, :3], Tc[1, 3:6], Tc[2, 6:] = -1., -1., -1.
    out = (np.zeros((3, F.DOF)), Tc, np.array([sd['V'], sd['Om'], sd['nu']]), 0, np.inf, E)
    if S.MULT_MODE == 'free':
        import scipy.sparse as sp
        out = out + (sp.csr_matrix((3, 0)),)
    return out


def franka_rows(q, links, obst, state, eta=F.ACT):
    """The travel bounds behind the skip threshold and the domain collar hold only under the caps, so every
    link gets its caps at every step, also when none of its boxes is active."""
    T, Z, O = F.fk(q)
    res, caps = [], {}
    for L in links:
        R, p = T[L['frame']][:3, :3], T[L['frame']][:3, 3]
        EL = _link_twist_map(q, L, (T, Z, O))
        sd = caps_for(q, L, state.u, EL)
        caps[L['frame']] = (sd['V'], sd['Om'], sd['nu'], sd['travel'])
        eta_L = max(eta, CONST['G'] * sd['travel'])
        capped = False
        for o in obst.values():
            if np.any(p + L['rad'] < o['dlo']) or np.any(p - L['rad'] > o['dhi']):
                continue
            r = S3.pair_rows(L, R, p, o, S3.I3, np.zeros(3), [(Z, O, L['n_joints'], 1., 0)], eta_L,
                             S3.UMAX7, sd=sd)
            capped |= r is not None             # pair rows carry the caps of L
            res.append(r)
        if not capped:
            res.append(_cap_block(EL, sd))
    state.caps = caps
    out = S3.stack(res, F.DOF)
    if len(out) == 6:                    # no active pair: keep the 'free' signature
        import scipy.sparse as sp
        out = out + (sp.csr_matrix((0, 0)),)
    return out


def solve_or_stop(u_nom, A, T, C, E, G_in, h_in, z_prev, Wc):
    """Hard QP; if it is infeasible, apply u = 0 (driftless: every pose, hence every h, stays constant)."""
    if Wc is not None and Wc.shape[1] == 0:
        Wc = None
    u, s, z = S.solve(u_nom, A, T, C, E, G_in, h_in, z_prev, Wc, soft=False)
    if s > 0:
        return np.zeros_like(u_nom), False, None
    return u, True, z


def describe():
    return dict(grow=GROW, floors=[FLOOR_V, FLOOR_OM, FLOOR_NU], speed_max_m_s=SPEED_MAX,
                G_max=CONST.get('G'), collar_bound_m=CONST.get('collar'), travel_max_m=CONST.get('travel_max'),
                lever_bounds_m={str(k): v for k, v in CONST.get('D', {}).items()})
