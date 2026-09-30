"""
Swarm escaping a maze: dynamic team B-spline SDF vs one static maze B-spline SDF,
kept apart by Bernstein-coefficient CBF constraints (same machinery as the 10 vs 10 script).

    Team shape:   S_T = {phi_T <= l_T}   (union of the robots' polygons, refitted every step)
    Maze shape:   S_M = {phi_M <= l_M}   (union of all wall rectangles, fitted once)
    Barrier:      g(x) = phi_T(x) + phi_M(x) - l_T - l_M - m  >= 0  for all x
                  => S_T and S_M never intersect.
The maze side is static (W_M_dot = 0, l_M constant), so only the team's terms move.

Nominal input: gradient descent on a geodesic (Dijkstra) distance-to-goal field computed once
on a grid of the maze inflated by the robot size; the CBF-QP alone is responsible for safety.
"""
import argparse
import time

import numpy as np
import torch
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

from bspline_sdf_bernstein_cbf_10v10 import (
    SHAPES, BSplineSDF2D, RobotShapes, adaptive_level, bernstein_shape_rows, intra_team_rows,
    robot_gap, rotate, solve_ldp_qp, speed_polygon)

def make_walls(domain=7.5, door=1.5, baffle_gap=1.2):
    """
    Maze walls as axis-aligned rectangles (x0, x1, y0, y1), thickness 0.2 m, in a 10 m box.
    Entrance / exit: doors of width `door` in the left wall centred at y = +-3.75.
    Four serpentine corridors, each with a baffle leaving an opening of `baffle_gap`.
    Fins at y = 0 outside the maze (to the domain edge) stop the team from walking around it.
    """
    e0, e1 = 3.75 - door / 2, 3.75 + door / 2
    return [
        (-5.1, 5.1, 4.9, 5.1),                      # outer top
        (-5.1, 5.1, -5.1, -4.9),                    # outer bottom
        (4.9, 5.1, -4.9, 4.9),                      # outer right
        (-5.1, -4.9, e1, 4.9),                      # outer left, above the entrance
        (-5.1, -4.9, -e0, e0),                      # outer left, between entrance and exit
        (-5.1, -4.9, -4.9, -e1),                    # outer left, below the exit
        (-4.9, 2.8, 2.3, 2.5),                      # H1, gap on the right
        (-2.8, 4.9, -0.1, 0.1),                     # H2, gap on the left
        (-4.9, 2.8, -2.5, -2.3),                    # H3, gap on the right
        (-1.1, -0.9, 2.5 + baffle_gap, 4.9),        # baffle, top corridor
        (0.9, 1.1, 0.1, 2.3 - baffle_gap),          # baffle, second corridor
        (-1.1, -0.9, -2.3 + baffle_gap, -0.1),      # baffle, third corridor
        (0.9, 1.1, -4.9, -2.5 - baffle_gap),        # baffle, bottom corridor
        (-domain, -5.1, -0.1, 0.1),                 # fin left
        (5.1, domain, -0.1, 0.1),                   # fin right
    ]


def wall_shapes(WALLS, device):
    polys, centers = [], []
    for x0, x1, y0, y1 in WALLS:
        c = np.array([(x0 + x1) / 2, (y0 + y1) / 2])
        polys.append(np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]) - c)
        centers.append(c)
    return RobotShapes(polys, device), np.array(centers)


def geodesic_fields(walls, wall_pos, goals, lo, hi, res, clearance):
    """
    Dijkstra distance-to-goal on an 8-connected grid of cells whose true maze SDF exceeds
    `clearance`, one field per goal. Cells inside the inflated walls take the value of their
    nearest free cell plus the distance to it, so the gradient points out of the walls too.
    Returns (x1d, fields (n_goals, res, res) indexed [ix, iy]).
    """
    x1d = np.linspace(lo, hi, res)
    gx, gy = np.meshgrid(x1d, x1d, indexing='ij')
    pts = torch.as_tensor(np.stack([gx.ravel(), gy.ravel()], axis=1), device=walls.device)
    sdf = walls.team_sdf(pts, wall_pos)[0].cpu().numpy().reshape(res, res)
    free = sdf > clearance

    idx = -np.ones((res, res), dtype=int)
    idx[free] = np.arange(free.sum())
    rows, cols, w = [], [], []
    for dx, dy in [(1, 0), (0, 1), (1, 1), (1, -1)]:
        a = idx[max(0, -dx):res - max(0, dx), max(0, -dy):res - max(0, dy)]
        b = idx[max(0, dx):res - max(0, -dx) or None, max(0, dy):res - max(0, -dy) or None]
        ok = (a >= 0) & (b >= 0)
        rows.append(a[ok])
        cols.append(b[ok])
        w.append(np.full(ok.sum(), np.hypot(dx, dy) * (x1d[1] - x1d[0])))
    n = free.sum()
    graph = coo_matrix((np.concatenate(w), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))

    src = [idx[np.argmin(np.abs(x1d - g[0])), np.argmin(np.abs(x1d - g[1]))] for g in goals]
    assert min(src) >= 0, "a goal lies inside the inflated walls"
    D_free = dijkstra(graph, directed=False, indices=src)                    # (n_goals, n)

    dist, (ix, iy) = ndimage.distance_transform_edt(~free, return_indices=True)
    fields = np.empty((len(goals), res, res))
    for k in range(len(goals)):
        D = np.full((res, res), np.inf)
        D[free] = D_free[k]
        fields[k] = D[ix, iy] + dist * (x1d[1] - x1d[0])
    return x1d, fields


def nominal_input(pos, goals, x1d, fields, u_max, k_goal, r_switch=0.6):
    """Follow -grad of each robot's geodesic field; plain goal attraction near the goal."""
    h = x1d[1] - x1d[0]
    u = np.zeros_like(pos)
    for k, (p, g) in enumerate(zip(pos, goals)):
        d = g - p
        if np.linalg.norm(d) < r_switch:
            u[k] = k_goal * d
        else:
            fx = np.clip((p[0] - x1d[0]) / h, 1, len(x1d) - 2)
            fy = np.clip((p[1] - x1d[0]) / h, 1, len(x1d) - 2)
            ix, iy = int(round(fx)), int(round(fy))
            F = fields[k]
            grad = np.array([F[ix + 1, iy] - F[ix - 1, iy], F[ix, iy + 1] - F[ix, iy - 1]])
            u[k] = -grad / max(np.linalg.norm(grad), 1e-9) * u_max
        s = np.linalg.norm(u[k])
        if s > u_max:
            u[k] *= u_max / s
    return u


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--frames', type=int, default=3000)
    ap.add_argument('--headless', action='store_true')
    ap.add_argument('--gif', action='store_true', help='save a GIF (every --gif-every frames)')
    ap.add_argument('--gif-every', type=int, default=5)
    ap.add_argument('--plot-every', type=int, default=5)
    ap.add_argument('--scale', type=float, default=1.0, help='robot shape size factor')
    ap.add_argument('--out', default='bspline_bernstein_maze.gif')
    ap.add_argument('--domain', type=float, default=7.5, help='half-width of the square domain [m]')
    ap.add_argument('--door', type=float, default=1.5, help='entrance / exit width [m]')
    ap.add_argument('--baffle-gap', type=float, default=1.2, help='opening left by each baffle [m]')
    args = ap.parse_args()

    import matplotlib
    if args.headless:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    domain = args.domain
    model = BSplineSDF2D(-domain, domain, n_cells=int(round(12 * domain)),        # h = 1/6 m
                         n_samples=int(round(320 * domain / 7.5)), lam=1e-3)

    # ---------- static maze: one B-spline SDF, fitted once ----------
    WALLS = make_walls(domain, args.door, args.baffle_gap)
    walls, wall_pos = wall_shapes(WALLS, model.dev)
    W_M, _ = model.fit_team(walls, wall_pos)
    pts_M = walls.samples + wall_pos[walls.samples_owner]
    l_M = model.predict(pts_M, W_M).max()                     # smallest enclosing level
    dW_M = torch.zeros(0, 2, model.n_ctrl, model.n_ctrl, device=model.dev, dtype=torch.float64)
    dv_M = np.zeros((1, 0, 2))

    # ---------- swarm ----------
    names = ['star', 'L', 'square', 'cross', 'triangle', 'U', 'hexagon', 'star', 'L', 'cross']
    rng = np.random.default_rng(1)
    team = RobotShapes([rotate(args.scale * SHAPES[s], rng.uniform(0, 2 * np.pi)) for s in names],
                       model.dev)
    N = len(names)
    sp = 0.85 * args.scale
    x0 = -domain + 0.75
    pos = np.array([[x0 + sp * (k % 2), 3.75 + sp * (k // 2 - 2)] for k in range(N)])
    goals = np.array([[x0 + sp * (k % 2), -3.75 + sp * (k // 2 - 2)] for k in range(N)])

    dt = 0.02
    u_max = 1.5
    k_goal = 2.0
    gamma = 6.0
    shape_margin = 0.05
    gamma_in = 6.0
    margin_in = 0.02
    eps_t = 1e-4

    t0 = time.perf_counter()
    x_plan, fields = geodesic_fields(walls, wall_pos, goals, -domain, domain, 301,
                                     clearance=team.bound_r.max() + 0.12)
    print(f"cells {model.K}x{model.K}, h = {model.h:.3f} m | l_maze = {l_M:+.4f} | "
          f"geodesic fields {time.perf_counter() - t0:.1f} s | path length "
          f"{np.mean([fields[k][np.argmin(np.abs(x_plan - p[0])), np.argmin(np.abs(x_plan - p[1]))] for k, p in enumerate(pos)]):.1f} m")

    G_spd, h_spd = speed_polygon(N, u_max)
    pad = lambda M: np.hstack([M, np.zeros((M.shape[0], 2))])

    vis_res = 300
    x_vis = np.linspace(-domain, domain, vis_res)
    X, Y = np.meshgrid(x_vis, x_vis, indexing='ij')
    zM_fine = model.predict_grid(np.linspace(-domain, domain, 600), W_M) - l_M
    pM_vis = model.predict_grid(x_vis, W_M)                   # fitted maze SDF (static)

    if not args.headless:
        plt.ion()
    fig, ax = plt.subplots(figsize=(8, 8))
    gif_frames = []
    worst = dict(gap=np.inf, sep=np.inf, beta=np.inf, slack=0.0)
    ctrl_ms = []
    z_prev = None
    arrived_at = None

    for frame in range(args.frames):
        t_start = time.perf_counter()
        W_T, dW_T = model.fit_team(team, pos)
        l_T, dv_T = adaptive_level(model, W_T, dW_T, team, pos)
        G_b, h_b, G_e, h_e, beta, n1 = bernstein_shape_rows(
            model, W_T, dW_T, l_T, dv_T, W_M, dW_M, l_M, dv_M, shape_margin, gamma, u_max)
        G_in, h_in = intra_team_rows(pos, 0, N, team.bound_r, gamma_in, margin_in)

        u_nom = nominal_input(pos, goals, x_plan, fields, u_max, k_goal)
        G = np.vstack([G_b, pad(G_in), G_e, pad(G_spd)])
        G[:, -2:] /= np.sqrt(eps_t)
        h = np.concatenate([h_b, h_in, h_e, h_spd])
        z, slack = solve_ldp_qp(np.concatenate([u_nom.ravel(), [0.0, 0.0]]), G, h,
                                len(h_b) + len(h_in), z_prev)
        z_prev = z
        u = z[:2 * N].reshape(N, 2)
        pos_before = pos
        pos = pos + u * dt
        ctrl_ms.append(1e3 * (time.perf_counter() - t_start))

        # ---------- safety diagnostics (state before the step) ----------
        gap = robot_gap(team, pos_before, walls, wall_pos)
        zT = model.predict_grid(np.linspace(-domain, domain, 600), W_T) - l_T
        sep = np.maximum(zT, zM_fine).min()
        worst['gap'] = min(worst['gap'], gap)
        worst['sep'] = min(worst['sep'], sep)
        worst['beta'] = min(worst['beta'], beta.min())
        worst['slack'] = max(worst['slack'], slack)
        err = np.linalg.norm(pos - goals, axis=1).max()
        if arrived_at is None and err < 0.1:
            arrived_at = frame

        if frame % 100 == 0:
            print(f"frame {frame:5d} | rows {len(h_b):5d} | min beta {beta.min():+.3f} | "
                  f"shape sep {sep:+.3f} | robot-wall gap {gap:+.3f} | slack {slack:.1e} | "
                  f"l_T {l_T:+.3f} | goal err {err:.2f} | ctrl {ctrl_ms[-1]:.1f} ms")

        want_gif = args.gif and frame % args.gif_every == 0
        if (not args.headless and frame % args.plot_every == 0) or want_gif:
            ax.clear()
            for (x0, x1, y0, y1) in WALLS:
                ax.add_patch(Polygon([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], closed=True,
                                     color='dimgray'))
            ax.contourf(X, Y, pM_vis, levels=[pM_vis.min() - 1, l_M], colors=['gray'], alpha=0.35)
            ax.contour(X, Y, pM_vis, levels=[l_M], colors='black', linewidths=1.0)
            pT = model.predict_grid(x_vis, W_T)
            ax.contourf(X, Y, pT, levels=[pT.min() - 1, l_T], colors=['royalblue'], alpha=0.3)
            ax.contour(X, Y, pT, levels=[l_T], colors='blue', linewidths=1.2)
            for poly, p in zip(team.polys, pos):
                ax.add_patch(Polygon(poly + p, closed=True, color='midnightblue', alpha=0.9))
            ax.plot(goals[:, 0], goals[:, 1], 'x', color='blue', ms=6)
            ax.set_title(f"Maze: team B-spline SDF (blue) vs fitted maze B-spline SDF (black) | frame {frame} | "
                         f"shape gap {sep:+.3f}", fontsize=10)
            ax.set_xlim(-domain, domain)
            ax.set_ylim(-domain, domain)
            ax.set_aspect('equal')
            if not args.headless:
                plt.pause(0.001)
            if want_gif:
                from PIL import Image
                fig.canvas.draw()
                w, hh = fig.canvas.get_width_height()
                buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(hh, w, 4)
                gif_frames.append(Image.fromarray(buf, 'RGBA').convert('RGB'))

        if arrived_at is not None and frame >= arrived_at + 50:
            break

    if args.gif:
        gif_frames[0].save(args.out, save_all=True, append_images=gif_frames[1:],
                           duration=40, loop=0)
        print(f">>> GIF saved: {args.out} ({len(gif_frames)} frames)")

    print("\n=== summary ===")
    print(f"min robot-wall gap (exact)     : {worst['gap']:+.4f} m")
    print(f"min team/maze shape separation : {worst['sep']:+.4f}  (>0: shapes never touch)")
    print(f"min Bernstein coefficient beta : {worst['beta']:+.4f}")
    print(f"max slack used                 : {worst['slack']:.2e}")
    print(f"all robots at goal (<0.1 m)    : "
          f"{'frame %d (t = %.1f s)' % (arrived_at, arrived_at * dt) if arrived_at is not None else 'no'}")
    print(f"max distance to goal           : {np.linalg.norm(pos - goals, axis=1).max():.3f} m")
    print(f"ctrl time median / p95         : {np.median(ctrl_ms):.1f} / {np.percentile(ctrl_ms, 95):.1f} ms")
    if not args.headless:
        plt.ioff()
        plt.show()


if __name__ == '__main__':
    main()
