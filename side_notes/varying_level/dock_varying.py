"""
Side note (not in the paper): spatially varying level lambda(q) instead of a single certified level.

Docking with a single certified level l* vs a level lambda(q) built from the per-cell certified
levels l_c.  Same B-spline field, same controller; only the level differs.

    python side_notes/varying_level/dock_varying.py
        -> side_notes/varying_level/dock_varying.gif, dock_varying_final.png, result.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

import cspace_sdf_cbf_compare as C                                        # noqa: E402
from cspace_experiments import CACHE, DOCK2_GOAL, DOCK2_START, dock2_setup, simulate_team  # noqa: E402
from dock_shapes import MOUTH_X, NOSE_TIP_X, SEAT_DEPTH                  # noqa: E402

import matplotlib                                                        # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                          # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter             # noqa: E402
from matplotlib.patches import Polygon                                   # noqa: E402

STEPS, DT, EVERY, HOLD = 900, 0.01, 8, 25


FLOOR = 0.0     # per-cell levels are clipped from below (always safe: lambda >= clipped >= l_c).
#                 Without clipping, deep-interior cells (l_c down to -1.2 m) next to boundary cells
#                 (+10 mm) force huge corrections and lambda becomes far larger than l* (tried:
#                 nose depth 0.139 m instead of 0.509 m).  Floors 0 / -10 / -30 mm gave depths
#                 0.519 / 0.502 / 0.480 m.


def varying_level(field, lc):
    """
    Spatially varying level lambda(q) in the field's own basis with lambda >= l_c on every cell c
    met by M (lc[c] = -inf elsewhere; those cells get the smallest finite l_c).  Control values
    start from l_c spread over the cells in which a control is interior; one correction pass then
    raises all controls of every cell whose Bernstein minimum is still below l_c by the deficit
    (raising only increases the Bernstein coefficients of the other cells).  Returns the control
    array Lam (same shape as field.W) and the min over cells of (Bernstein min of lambda - l_c).
    """
    L = np.where(np.isfinite(lc), lc, lc[np.isfinite(lc)].min())
    bases = (field.bx, field.by, field.bt)

    def spread(V, offsets):
        for ax, b in enumerate(bases):
            idx = b.ctrl_idx(np.arange(b.K))                    # (K, 4) controls of every cell
            out = np.full(V.shape[:ax] + (b.n,) + V.shape[ax + 1:], -np.inf)
            Vm = np.moveaxis(V, ax, 0)
            om = np.moveaxis(out, ax, 0)
            for o in offsets:
                np.maximum.at(om, idx[:, o], Vm)
            V = np.moveaxis(om, 0, ax)
        return V

    def bern_min(Lam):
        f = C.Field3D(field.kind, field.bx.hi, field.bx.K, field.bt.K)
        f.set_W(Lam)
        return f.C.min(axis=(3, 4, 5))

    Lam = spread(L, (1, 2))
    Lam = np.where(np.isfinite(Lam), Lam, L.min())
    deficit = np.maximum(L - bern_min(Lam), 0.0)
    Lam = Lam + np.maximum(spread(deficit, (0, 1, 2, 3)), 0.0)
    return Lam, float((bern_min(Lam) - L).min())


def varying_field(A, B, f, D):
    """psi = phi - lambda (same basis) with its certified level; cached."""
    fn = os.path.join(CACHE, f'dock_craft_varying_floor{FLOOR * 1e3:+.0f}mm.npz')
    if os.path.exists(fn):
        d = np.load(fn)
        fv = C.Field3D('bspline', D, f.bx.K, f.bt.K)
        fv.set_W(d['W'])
        return fv, float(d['l']), json.loads(str(d['meta']))
    cert = C.Certifier(f, A, B)
    fn_lc = os.path.join(CACHE, 'dock_craft_lc.npy')
    if os.path.exists(fn_lc):
        lc, ok_c, t_c = np.load(fn_lc), None, None
    else:
        lc, ok_c, t_c = cert.certify_cells(tol=2e-3)
        np.save(fn_lc, lc)
    Lam, margin = varying_level(f, np.maximum(lc, FLOOR))
    fv = C.Field3D('bspline', D, f.bx.K, f.bt.K)
    fv.set_W(f.W - Lam)
    cv = C.Certifier(fv, A, B)
    l_v, ok_v, _, t_v = cv.certify(verbose=False)
    band_v, band_min_v = cv.certify_band(l_v)
    reg_v, _, _ = cv.certify_regular(l_v)
    fin = lc[np.isfinite(lc)]
    meta = dict(floor=FLOOR, lam_min=float(Lam.min()), lam_max=float(Lam.max()),
                cells_met=int(fin.size), lc_min=float(fin.min()), lc_med=float(np.median(fin)),
                lc_max=float(fin.max()), cells_certified=ok_c, t_cells=t_c, bern_margin=margin,
                l_star=l_v, certified=ok_v, t_cert=t_v, band=band_v, band_min=band_min_v, regular=reg_v)
    np.savez(fn, W=fv.W, l=l_v, meta=json.dumps(meta))
    return fv, l_v, meta


def run(GA, GB, A, B, field):
    L = simulate_team([GA, GB], [A, B], {(0, 1): field}, DOCK2_START, DOCK2_GOAL, 'cspace',
                      steps=STEPS, v_max=np.array([0.0, 1.0]), w_max=np.array([0.0, 2.0]), record=True)
    traj = np.array(L['traj'])[:, 1]
    tip = traj[:, :2] + NOSE_TIP_X * np.stack([np.cos(traj[:, 2]), np.sin(traj[:, 2])], axis=1)
    gap = np.array([C.true_distance(GA, GB @ C.rot(x[2]).T + x[:2]) for x in traj])
    return traj, MOUTH_X - tip[:, 0], gap, L


def main():
    GA, GB, A, B, fields, meta = dock2_setup()
    f, D, l_star = fields['bspline']
    fv, l_v, meta_v = varying_field(A, B, f, D)
    print('varying level:', meta_v)
    runs = [('single level $l^\\star$ = %.1f mm' % (1e3 * l_star), '#1f77b4',
             run(GA, GB, A, B, (f, D, l_star))),
            ('per-patch level $\\lambda(q)$', '#d62728', run(GA, GB, A, B, (fv, D, l_v)))]
    res = dict(varying=meta_v, l_star=l_star)
    for key, (name, _, (traj, depth, gap, L)) in zip(('single', 'varying'), runs):
        res[key] = dict(depth=float(depth[-1]), shortfall=float(SEAT_DEPTH - depth[-1]),
                        final_gap=float(gap[-1]), min_gap=float(gap.min()),
                        min_h=float(np.nanmin(L['h_t'])), tv=float(L['tv']))
        print(f"{name}: {res[key]}")
    json.dump(res, open(os.path.join(HERE, 'result.json'), 'w'), indent=2)

    plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
    fig = plt.figure(figsize=(10, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.4], hspace=0.28, wspace=0.06)
    texts, polys_gt, polys_fit, trails = [], [], [], []
    for k, (name, col, (traj, _, _, _)) in enumerate(runs):
        ax = fig.add_subplot(gs[0, k])
        ax.add_patch(Polygon(GA, closed=True, color='0.6', lw=0))
        ax.add_patch(Polygon(A.dense, closed=True, fill=False, ec='k', lw=0.9))
        pg = ax.add_patch(Polygon(GB, closed=True, color='0.6', alpha=0.8, lw=0))
        pf = ax.add_patch(Polygon(B.dense, closed=True, fill=False, ec=col, lw=1.2))
        tr, = ax.plot([], [], color=col, lw=0.8)
        ax.set_xlim(-1.0, 4.3)
        ax.set_ylim(-1.1, 1.25)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(name, color=col)
        az = ax.inset_axes([0.60, 0.02, 0.39, 0.36])                  # zoom on the receptacle
        az.add_patch(Polygon(GA, closed=True, color='0.6', lw=0))
        az.add_patch(Polygon(A.dense, closed=True, fill=False, ec='k', lw=0.7))
        zg = az.add_patch(Polygon(GB, closed=True, color='0.6', alpha=0.8, lw=0))
        zf = az.add_patch(Polygon(B.dense, closed=True, fill=False, ec=col, lw=1.0))
        az.set_xlim(MOUTH_X - SEAT_DEPTH - 0.08, MOUTH_X + 0.12)
        az.set_ylim(-0.2, 0.2)
        az.set_aspect('equal')
        az.set_xticks([])
        az.set_yticks([])
        for s in az.spines.values():
            s.set_color(col)
        texts.append(ax.text(0.02, 0.97, '', transform=ax.transAxes, va='top', fontsize=8.5))
        polys_gt.append((pg, zg))
        polys_fit.append((pf, zf))
        trails.append(tr)
    axd = fig.add_subplot(gs[1, :])
    t = np.arange(STEPS) * DT
    lines = []
    for name, col, (_, depth, _, _) in runs:
        axd.plot(t[:len(depth)], depth, color=col, lw=0.6, alpha=0.25)
        ln, = axd.plot([], [], color=col, lw=1.6, label=name)
        lines.append(ln)
    axd.axhline(0, color='k', lw=0.5, ls=':')
    axd.axhline(SEAT_DEPTH, color='k', lw=0.5, ls='--')
    axd.text(t[-1], SEAT_DEPTH + 0.03, 'full seating', ha='right', fontsize=8)
    axd.text(t[-1], 0.03, 'mouth', ha='right', fontsize=8)
    axd.set_xlim(0, t[-1])
    axd.set_ylim(-2.5, SEAT_DEPTH + 0.2)
    axd.set_xlabel('$t$ [s]')
    axd.set_ylabel('nose depth past the mouth [m]')
    axd.legend(loc='center right', frameon=False)

    n_frames = STEPS // EVERY + HOLD

    def update(fr):
        j = min(fr * EVERY, STEPS - 1)
        for k, (name, col, (traj, depth, gap, _)) in enumerate(runs):
            jj = min(j, len(traj) - 1)
            x = traj[jj]
            for p in polys_gt[k]:
                p.set_xy(GB @ C.rot(x[2]).T + x[:2])
            for p in polys_fit[k]:
                p.set_xy(B.world(x))
            trails[k].set_data(traj[:jj + 1, 0], traj[:jj + 1, 1])
            texts[k].set_text(f"nose depth {depth[jj]:+.3f} m / {SEAT_DEPTH:.2f} m\n"
                              f"true gap {1e3 * gap[jj]:5.1f} mm")
            lines[k].set_data(t[:jj + 1], depth[:jj + 1])
        return []

    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    anim.save(os.path.join(HERE, 'dock_varying.gif'), writer=PillowWriter(fps=20), dpi=80)
    update(n_frames - 1)
    fig.savefig(os.path.join(HERE, 'dock_varying_final.png'), dpi=110, bbox_inches='tight')
    print('>>> side_notes/varying_level/dock_varying.gif')


if __name__ == '__main__':
    main()
