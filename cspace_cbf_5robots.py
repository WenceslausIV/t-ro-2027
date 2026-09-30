"""
Five robots with random non-convex shapes swap positions on a circle,
using pairwise C-space B-spline SDF CBFs with certified enclosures. No deadlock resolution.

    python cspace_cbf_5robots.py                 # five random shapes
    python make_paper_figs.py five               # regenerate the corresponding paper figure

Pipeline (all offline except the last line):
  1. ground-truth shape  G_i : random radial Fourier polygon
  2. fitted curve        p_i : periodic cubic B-spline fitted to outward-offset samples
                               (ICRA'27 style); certified  G_i  inside  p_i  by Bezier control-
                               hull tests with de Casteljau subdivision + winding
  3. for every pair (i, j): C-space SDF phi_ij(t, theta) of the fitted shapes (3-D B-spline),
                               certified level l*_ij >= max over {b(tau) - R(theta)a(sigma)}
  4. online: h_ij = phi_ij(q_ij) - l*_ij >= 0 for all pairs, one joint CBF-QP.
Safety chain: h_ij >= 0  =>  fitted shapes never touch  =>  ground-truth shapes never touch.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path as FilePath
import time

import numpy as np
from matplotlib.path import Path

from cspace_sdf_cbf_compare import (ACT_DELTA, MB, MB_INV, TWO_PI, Basis1D, Certifier, Field3D, SmoothShape, bern3,
                                    cspace_sdf_samples, cspace_sdf_slices, rot, sample_slices, true_distance, wrap)
from sdf_cbf_utils import solve_ldp_qp


# ---------------------------------------------------------------------------------------------
# 1. random ground-truth shapes
# ---------------------------------------------------------------------------------------------
def random_shape(rng, n=400):
    """Star-shaped smooth non-convex polygon (CCW): r(phi) = r0 (1 + sum_k a_k cos(k phi + p_k))."""
    phi = np.arange(n) * TWO_PI / n
    while True:
        r0 = rng.uniform(0.3, 0.38)
        ks = np.arange(2, 6)
        a = rng.uniform(0.05, 0.28, len(ks)) * rng.choice([0, 1], len(ks), p=[0.3, 0.7])
        p = rng.uniform(0, TWO_PI, len(ks))
        r = r0 * (1 + (a[:, None] * np.cos(ks[:, None] * phi + p[:, None])).sum(0))
        P = np.stack([r * np.cos(phi), r * np.sin(phi)], axis=1)
        hull_area = convex_hull_area(P)
        area = 0.5 * np.abs(np.dot(P[:, 0], np.roll(P[:, 1], -1)) - np.dot(P[:, 1], np.roll(P[:, 0], -1)))
        if r.min() > 0.4 * r0 and area / hull_area < 0.9:                 # clearly non-convex
            return P


def convex_hull_area(P):
    pts = sorted(map(tuple, P))

    def half(seq):
        h = []
        for q in seq:
            while len(h) >= 2 and ((h[-1][0] - h[-2][0]) * (q[1] - h[-2][1])
                                   - (h[-1][1] - h[-2][1]) * (q[0] - h[-2][0])) <= 0:
                h.pop()
            h.append(q)
        return h
    H = np.array(half(pts)[:-1] + half(pts[::-1])[:-1])
    return 0.5 * np.abs(np.dot(H[:, 0], np.roll(H[:, 1], -1)) - np.dot(H[:, 1], np.roll(H[:, 0], -1)))


# ---------------------------------------------------------------------------------------------
# 2. fitted B-spline boundary with certified enclosure of the ground truth
# ---------------------------------------------------------------------------------------------
def fit_curve(G, K, offset, lam=1e-4, n_fit=200):
    idx = np.linspace(0, len(G), n_fit, endpoint=False).astype(int)
    P = G[idx]
    tang = np.roll(G, -1, 0)[idx] - np.roll(G, 1, 0)[idx]
    nrm = np.stack([tang[:, 1], -tang[:, 0]], axis=1)                    # outward for CCW
    P = P + offset * nrm / np.linalg.norm(nrm, axis=1, keepdims=True)
    seg = np.linalg.norm(np.diff(np.vstack([P, P[:1]]), axis=0), axis=1)
    tau = np.r_[0, np.cumsum(seg)[:-1]] / seg.sum()                      # chord-length parameter
    Phi = Basis1D('bspline', 0.0, 1.0, K, True).design(tau)
    return np.linalg.solve(Phi.T @ Phi + lam * np.eye(K), Phi.T @ P)     # control points (K, 2)


def _seg_cross(p1, p2, q1, q2):
    """Proper/touching intersection of segments p1p2 and q1q2 (broadcast)."""
    def orient(a, b, c):
        return (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])
    d1, d2 = orient(q1, q2, p1), orient(q1, q2, p2)
    d3, d4 = orient(p1, p2, q1), orient(p1, p2, q2)
    return (d1 * d2 <= 0) & (d3 * d4 <= 0)


def hulls_meet_polygon(H, G):
    """H (M,4,2) Bezier control points; True where conv(H[m]) intersects the polygon G."""
    path = Path(G)
    meet = path.contains_points(H.reshape(-1, 2)).reshape(-1, 4).any(1)  # hull vertex inside G
    tri = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]                   # conv(4 pts) = union
    for a, b, c in tri:
        A_, B_, C_ = H[:, a, None], H[:, b, None], H[:, c, None]
        def cr(u, v, w):
            return (v[..., 0] - u[..., 0]) * (w[..., 1] - u[..., 1]) - (v[..., 1] - u[..., 1]) * (w[..., 0] - u[..., 0])
        s1, s2, s3 = cr(A_, B_, G[None]), cr(B_, C_, G[None]), cr(C_, A_, G[None])
        inside = ((s1 >= 0) & (s2 >= 0) & (s3 >= 0)) | ((s1 <= 0) & (s2 <= 0) & (s3 <= 0))
        meet |= inside.any(1)                                            # G vertex inside hull
    G2 = np.roll(G, -1, 0)
    for a, b in itertools.combinations(range(4), 2):                     # hull edges vs G edges
        meet |= _seg_cross(H[:, a, None], H[:, b, None], G[None], G2[None]).any(1)
    return meet


def sub_matrix_np(s0, s1):
    """(R,4,4): cubic Bernstein coefficients on [0,1] -> on the sub-interval [s0, s1]."""
    d = s1 - s0
    T = np.zeros((len(s0), 4, 4))
    for j in range(4):
        for k in range(j, 4):
            T[:, j, k] = math.comb(k, j) * s0 ** (k - j) * d ** j
    return MB_INV @ T @ MB


def winding_number(poly, q):
    ang = np.arctan2(poly[:, 1] - q[1], poly[:, 0] - q[0])
    return int(round(wrap(np.diff(np.r_[ang, ang[:1]])).sum() / TWO_PI))


def certify_curve_encloses(ctrl, G, max_depth=6):
    """
    Lemma (ICRA'27): if the control hulls of a finite subdivision of the curve's spans all miss G
    and the polygon through the pieces' endpoints (in curve order) has nonzero winding around a
    point of G, then G lies inside the curve.
    """
    shape = SmoothShape(ctrl)
    k = np.arange(len(shape.bez))
    s0, s1 = np.zeros(len(k)), np.ones(len(k))
    acc_k, acc_s = [], []
    for depth in range(max_depth + 1):
        H = sub_matrix_np(s0, s1) @ shape.bez[k]
        bad = hulls_meet_polygon(H, G)
        acc_k.append(k[~bad])
        acc_s.append(s0[~bad])
        if not bad.any():
            break
        if depth == max_depth:
            return shape, False
        k, s0, s1 = k[bad], s0[bad], s1[bad]
        m = (s0 + s1) / 2
        k, s0, s1 = np.r_[k, k], np.r_[s0, m], np.r_[m, s1]
    kk, ss = np.concatenate(acc_k), np.concatenate(acc_s)
    order = np.lexsort((ss, kk))                                          # curve order
    chord_poly = np.einsum('ni,nid->nd', bern3(ss[order]), shape.bez[kk[order]])
    return shape, winding_number(chord_poly, G[0]) != 0


def enclose(G, K=24):
    for off in np.arange(0.004, 0.08, 0.004):
        shape, ok = certify_curve_encloses(fit_curve(G, K, off, n_fit=max(200, 4 * K)), G)
        if ok:
            return shape, off
    raise RuntimeError("no enclosing fit found")


# ---------------------------------------------------------------------------------------------
# 3. pairwise C-space SDFs
# ---------------------------------------------------------------------------------------------
def pair_field(A, B, n_th=96, n_xy=100):
    D = round(A.rho + B.rho + 0.3, 2)
    th = np.arange(n_th) * TWO_PI / n_th
    gxy = np.linspace(-D, D, n_xy)
    f = Field3D('bspline', D, int(round(2 * D / 0.0625)), 48)       # h ~ 6 cm
    f.fit(cspace_sdf_samples(A, B, th, D, 0.01, gxy, gxy), gxy, gxy, th)   # = slices + sampling, on the GPU
    cert = Certifier(f, A, B)
    l_star, ok, _, _ = cert.certify(verbose=False)
    band_ok, band_min = cert.certify_band(l_star)
    reg_ok, _, reg_open = cert.certify_regular(l_star)
    return f, D, l_star, ok and band_ok and reg_ok, (ok, band_ok, band_min, reg_ok, reg_open)


# ---------------------------------------------------------------------------------------------
def five_robot_setup(seed=3, scale=1.0, radius=2.6, jitter_seed=None, jitter=(0.0, 0.0), goal_shift=None,
                     n_robots=5):
    """Shared GIF/paper setup, with a geometry-specific cache separate from random statistics.
    scale multiplies the ground-truth shapes; with jitter_seed, start and goal angles on the circle
    of the given radius are perturbed by up to jitter = (start, goal) radians (breaks the symmetric
    deadlock of exact antipodal swaps)."""
    rng = np.random.default_rng(seed)
    GT = [scale * random_shape(rng) for _ in range(5)]
    counts = [24] * 5
    fits = [enclose(g, K=k) for g, k in zip(GT, counts)]
    shapes, offsets = [s for s, _ in fits], [o for _, o in fits]
    for i, (s, off) in enumerate(fits):
        print(f'robot {i}: {counts[i]} controls, enclosure passed, offset {off:.3f} m, rho {s.rho:.3f}', flush=True)
    root = FilePath(__file__).resolve().parent
    digest = hashlib.sha256(b'five-random-v1;96;100;0.0625;48;0.002')
    for s in shapes:
        digest.update(s.ctrl.tobytes())
    cache = root / 'cache' / f'five_random_{digest.hexdigest()[:16]}.npz'
    cache.parent.mkdir(exist_ok=True)
    # Reuse the original statistical cache only when both geometry and fits match exactly.
    legacy, legacy_meta = root / 'cache' / 'shapes5.npz', root / 'cache' / 'shapes5.json'
    if not cache.exists() and legacy.exists() and legacy_meta.exists():
        with np.load(legacy, allow_pickle=False) as saved:
            original = json.loads(legacy_meta.read_text(encoding='utf-8'))
            matches = (np.array_equal(saved['GT'], np.array(GT))
                       and np.array_equal(saved['ctrl'], np.array([s.ctrl for s in shapes]))
                       and len(original['pairs']) == 10
                       and all(m['ok'] for m in original['pairs'].values()))
            if matches:
                np.savez(cache, pairs=json.dumps(original['pairs']),
                         **{f'W_{key}': saved[f'W_{key}'] for key in original['pairs']})
    fields, pairs, data = {}, {}, {}
    if cache.exists():
        with np.load(cache, allow_pickle=False) as saved:
            pairs = json.loads(str(saved['pairs']))
            if not all(m['ok'] for m in pairs.values()):
                raise RuntimeError(f'Uncertified pair in {cache}')
            for key, m in pairs.items():
                f = Field3D('bspline', m['D'], m['Kxy'], m['Kth'])
                f.set_W(saved[f'W_{key}'])
                fields[tuple(map(int, key.split('_')))] = (f, m['D'], m['l'])
        print(f'loaded {cache.name}', flush=True)
    else:
        for i, j in itertools.combinations(range(5), 2):
            t0 = time.perf_counter()
            f, D, level, ok, info = pair_field(shapes[i], shapes[j])
            print(f'pair ({i},{j}): l* {level:+.4f}, enclosure/band/regularity {info}, '
                  f'{time.perf_counter()-t0:.1f} s', flush=True)
            if not ok:
                raise RuntimeError(f'Pair ({i},{j}) failed certification; refusing simulation')
            key = f'{i}_{j}'
            fields[i, j] = (f, D, level)
            pairs[key] = dict(D=D, Kxy=f.bx.K, Kth=f.bt.K, l=level, ok=bool(ok), band_min=info[2])
            data[f'W_{key}'] = f.W
        np.savez(cache, pairs=json.dumps(pairs), **data)
    ang = np.pi / 2 + TWO_PI * np.arange(5) / 5
    a_s, a_g = ang.copy(), (ang + np.pi if goal_shift is None else ang + TWO_PI * goal_shift / 5)
    if jitter_seed is not None:
        jr = np.random.default_rng(jitter_seed)
        a_s = a_s + jr.uniform(-jitter[0], jitter[0], 5)
        a_g = a_g + jr.uniform(-jitter[1], jitter[1], 5)
    # optional third jitter entry: start distances from the center vary by up to +/- jitter[2]
    r_s = radius + (jr.uniform(-jitter[2], jitter[2], 5) if jitter_seed is not None and len(jitter) > 2 else 0.)
    starts = np.stack([r_s*np.cos(a_s), r_s*np.sin(a_s), rng.uniform(-np.pi, np.pi, 5)], axis=1)
    gx, gy = ((radius*np.cos(a_g), radius*np.sin(a_g)) if (jitter_seed is not None or goal_shift is not None)
              else (-radius*np.cos(ang), -radius*np.sin(ang)))
    goals = np.stack([gx, gy, rng.uniform(-np.pi, np.pi, 5)], axis=1)
    if n_robots != 5:
        # first n_robots shapes (their pair fields are a subset of the cached ones); antipodal swap on
        # a circle, start/goal angles jittered by jitter[0]/jitter[1] and start/goal radii by jitter[2]
        n = n_robots
        j = list(jitter) + [0.] * (3 - len(jitter))
        jr = np.random.default_rng(jitter_seed if jitter_seed is not None else 0)
        on = 1. if jitter_seed is not None else 0.
        ang = np.pi / 2 + TWO_PI * np.arange(n) / n
        a_s = ang + on * jr.uniform(-j[0], j[0], n)
        a_g = ang + np.pi + on * jr.uniform(-j[1], j[1], n)
        r_s = radius + on * jr.uniform(-j[2], j[2], n)
        r_g = radius + on * jr.uniform(-j[2], j[2], n)
        hd = np.random.default_rng(seed + 100).uniform(-np.pi, np.pi, (2, n))
        starts = np.stack([r_s * np.cos(a_s), r_s * np.sin(a_s), hd[0]], axis=1)
        goals = np.stack([r_g * np.cos(a_g), r_g * np.sin(a_g), hd[1]], axis=1)
        GT, shapes, offsets = GT[:n], shapes[:n], offsets[:n]
        fields = {k: v for k, v in fields.items() if max(k) < n}
    meta = dict(shape_set='random', seed=seed, offsets=offsets, pairs=pairs, cache=cache.name,
                names=[f'random {i}' for i in range(5)])
    return GT, shapes, fields, starts, goals, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--steps', type=int, default=1500)
    ap.add_argument('--dt', type=float, default=0.01)
    ap.add_argument('--gif-every', type=int, default=6)
    ap.add_argument('--out', default='cspace_cbf_5robots.gif')
    args = ap.parse_args()
    N = 5
    GT, shapes, fields, x, goals, meta = five_robot_setup(args.seed)

    dt, v_max, w_max, gamma, k_v, k_w = args.dt, 1.0, 2.0, 5.0, 1.5, 2.0
    nv = 3 * N
    ang12 = np.linspace(0, TWO_PI, 12, endpoint=False)
    nk = np.stack([np.cos(ang12), np.sin(ang12)], axis=1)
    G_box, h_box = [], []
    for i in range(N):
        for n_ in nk:
            r = np.zeros(nv)
            r[3 * i:3 * i + 2] = -n_
            G_box.append(r)
            h_box.append(-v_max * np.cos(np.pi / 12))
        for sgn in (1, -1):
            r = np.zeros(nv)
            r[3 * i + 2] = -sgn
            G_box.append(r)
            h_box.append(-w_max)
    G_box, h_box = np.array(G_box), np.array(h_box)
    J2 = np.array([[0.0, -1.0], [1.0, 0.0]])

    log = dict(x=[], h=[], d=[], t=[])
    z_prev = None
    for k in range(args.steps):
        t0 = time.perf_counter()
        u_nom = np.zeros(nv)
        for i in range(N):
            v = k_v * (goals[i, :2] - x[i, :2])
            s = np.linalg.norm(v)
            u_nom[3 * i:3 * i + 2] = v if s <= v_max else v * v_max / s
            u_nom[3 * i + 2] = np.clip(k_w * wrap(goals[i, 2] - x[i, 2]), -w_max, w_max)
        rows, rhs, hs = [], [], []
        for (i, j), (f, D, l_star) in fields.items():
            t = rot(-x[j, 2]) @ (x[i, :2] - x[j, :2])
            if np.any(np.abs(t) >= D - ACT_DELTA):
                continue
            th = (x[i, 2] - x[j, 2]) % TWO_PI
            phi, g = f.eval(np.r_[t, th], grad=True)
            h = phi[0] - l_star
            gt, gth = g[0, :2], g[0, 2]
            r = np.zeros(nv)
            r[3 * i:3 * i + 2] = rot(x[j, 2]) @ gt
            r[3 * j:3 * j + 2] = -rot(x[j, 2]) @ gt
            r[3 * i + 2] = gth
            r[3 * j + 2] = -gt @ (J2 @ t) - gth
            rows.append(r)
            rhs.append(-gamma * h)
            hs.append(h)
        G = np.vstack(rows + [G_box]) if rows else G_box
        hh = np.r_[rhs, h_box]
        u, _ = solve_ldp_qp(u_nom, G, hh, len(rows), z_prev)
        z_prev = u
        log['t'].append(time.perf_counter() - t0)

        world = [GT[i] @ rot(x[i, 2]).T + x[i, :2] for i in range(N)]
        dmin = min(true_distance(world[i], world[j]) for i, j in fields)   # ground-truth shapes
        log['x'].append(x.copy())
        log['h'].append(min(hs) if hs else np.nan)
        log['d'].append(dmin)
        x = x + dt * u.reshape(N, 3)
    for key in log:
        log[key] = np.array(log[key])

    err = np.linalg.norm(log['x'][-1, :, :2] - goals[:, :2], axis=1)
    print("\n=== summary ===")
    print(f"min ground-truth distance over all pairs : {log['d'].min():+.4f} m")
    print(f"min h over all active pairs             : {np.nanmin(log['h']):+.4f}")
    print(f"robots within 5 cm of goal              : {(err < 0.05).sum()} / {N}   "
          f"(final errors {np.round(err, 2)})")
    print(f"online CBF+QP per step, median / max    : {1e3 * np.median(log['t']):.2f} / "
          f"{1e3 * log['t'].max():.2f} ms")

    result_dir = FilePath(__file__).resolve().parent / 'results'
    result_dir.mkdir(exist_ok=True)
    stem = 'five_random'
    np.savez(result_dir / f'{stem}_trajectory.npz', **log, final=x, goals=goals, dt=dt,
             **{f'GT_{i}': g for i, g in enumerate(GT)},
             **{f'ctrl_{i}': s.ctrl for i, s in enumerate(shapes)})
    meta.update(min_gt=float(log['d'].min()), min_h=float(np.nanmin(log['h'])),
                final_errors=err.tolist(), t_med_ms=float(1e3*np.median(log['t'])),
                t_max_ms=float(1e3*log['t'].max()), steps=args.steps, dt=dt)
    position_error = np.linalg.norm(log['x'][:, :, :2] - goals[None, :, :2], axis=2)
    angle_error = np.abs(wrap(log['x'][:, :, 2] - goals[None, :, 2]))
    reached = np.flatnonzero(((position_error < 0.05) & (angle_error < 0.05)).all(axis=1))
    meta['all_goals_first_time'] = float(reached[0] * dt) if len(reached) else None
    from polygon_validation import audit_trajectory
    meta['polygon_audit'] = audit_trajectory(GT, np.concatenate([log['x'], x[None]]))
    print(f"polygon intersection audit: {meta['polygon_audit']['intersection_count']} contacts "
          f"at {meta['polygon_audit']['checked_poses']} sampled poses", flush=True)
    (result_dir / f'{stem}.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')

    make_gif(log, GT, shapes, goals, args)


def make_gif(log, GT, shapes, goals, args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from PIL import Image

    N = len(GT)
    cols = plt.cm.tab10(np.arange(N))
    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1])
    ax, axp = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    n = len(log['d'])
    tt = np.arange(n) * args.dt
    frames = []
    for k in list(range(0, n, args.gif_every)) + [n - 1]:
        ax.clear()
        for i in range(N):
            ax.add_patch(Polygon(shapes[i].world(goals[i]), closed=True, fill=False, ec=cols[i],
                                 ls=':', lw=1))
            ax.plot(log['x'][:k + 1, i, 0], log['x'][:k + 1, i, 1], color=cols[i], lw=0.8, alpha=0.6)
            xi = log['x'][k, i]
            ax.add_patch(Polygon(GT[i] @ rot(xi[2]).T + xi[:2], closed=True, color='0.55'))
            ax.add_patch(Polygon(shapes[i].world(xi), closed=True, fill=False, ec=cols[i], lw=2))
        h = log['h'][k]
        ax.set_title(f"t = {k * args.dt:.2f} s | min h = {'--' if np.isnan(h) else f'{h:+.3f}'} | "
                     f"min true dist = {log['d'][k]:+.3f} m\n"
                     "gray: ground-truth shape, colored: certified fitted B-spline curve", fontsize=10)
        ax.set_xlim(-3.6, 3.6)
        ax.set_ylim(-3.6, 3.6)
        ax.set_aspect('equal')
        axp.clear()
        axp.plot(tt[:k + 1], log['d'][:k + 1], color='k', label='min ground-truth distance (all pairs)')
        axp.plot(tt[:k + 1], log['h'][:k + 1], color='tab:green', ls=':', label='min h (active pairs)')
        axp.axhline(0, color='0.3', lw=0.8)
        axp.set_xlim(0, tt[-1])
        axp.set_ylim(-0.05, 1.5)
        axp.set_xlabel('t [s]')
        axp.legend(fontsize=8, loc='upper right')
        # Fixed margins keep the two-line title inside the canvas on every frame.
        fig.subplots_adjust(left=0.055, right=0.985, bottom=0.09, top=0.90, wspace=0.19)
        fig.canvas.draw()
        w, hh = fig.canvas.get_width_height()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(hh, w, 4)
        frames.append(Image.fromarray(buf, 'RGBA').convert('RGB'))
    frames[0].save(args.out, save_all=True, append_images=frames[1:],
                   duration=max(10, round(args.dt * args.gif_every * 1000)), loop=0)
    print(f">>> saved {args.out} ({len(frames)} frames)")


if __name__ == '__main__':
    main()
