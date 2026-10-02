"""Side-by-side rotating 3D video of one star-tube goal that every listed method reaches.

    python prototype_3d/star_random_goals_video.py OUT.mp4 FOLDER[:LABEL] ... [--trial K] [--fps 20]
Each panel replays one method's saved trajectory on a common clock (a method that arrives early holds its final
pose); the camera turns once around the scene. Sphere methods also show their tube spheres (translucent).
Without --trial, the goal reached by all methods with the largest spread of arrival times is chosen.
"""
import argparse
import json

import numpy as np

import star_random_goals as G
import franka3d as F
from franka_gif import decimate, shade
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection         # noqa: E402
from star_random_goals_fig import sphere_mesh, ROBOT, TUBE, TSPH  # noqa: E402

GOAL = (.20, .70, .30)


def load(folder):
    recs = {}
    for p in sorted((G.OUT / folder).glob('trial_*.json')):
        m = json.loads(p.read_text())['metrics']
        recs[int(p.stem[-2:])] = m
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('methods', nargs='+')
    ap.add_argument('--trial', type=int, default=None)
    ap.add_argument('--fps', type=int, default=20)
    ap.add_argument('--hold', type=float, default=1.5, help='seconds to hold after the last arrival')
    a = ap.parse_args()
    spec = [(m.split(':')[0], m.split(':')[1] if ':' in m else m.split(':')[0]) for m in a.methods]
    recs = {f: load(f) for f, _ in spec}
    both = [k for k in range(30) if all(recs[f].get(k, {}).get('reached_s') is not None for f, _ in spec)]
    assert both, 'no goal is reached by every method'
    k = a.trial if a.trial is not None else max(both, key=lambda i: np.ptp([recs[f][i]['reached_s'] for f, _ in spec]))
    print('goals reached by all:', both, '-> trial', k, flush=True)

    q0, links, B, RB, pB = G.scene()
    Vt, Ft = G.SV.exact_surface(B)
    tube = (Vt @ RB.T + pB)[Ft]
    tube_col = np.c_[shade(tube, TUBE), np.full(len(tube), .7)]
    meshes = [decimate(L['V'], L['F'], .012) for L in links]
    qg = np.array(json.loads((G.OUT / 'goals.json').read_text())['goals'][k])
    traj = {f: np.load(G.OUT / f / f'trial_{k:02d}.npz')['q'] for f, _ in spec}
    spheres = {}
    for f, _ in spec:
        st = json.loads((G.OUT / f / 'setup.json').read_text())
        if st['method'] == 'spheres':
            Ct, rt = G.tube_spheres(B, st['voxel_mm'] / 1000, int(f.split('_T')[1]))
            spheres[f] = sphere_mesh(Ct @ RB.T + pB, rt, 8)
    T_end = max(recs[f][k]['reached_s'] for f, _ in spec) + a.hold
    n_frames = int(T_end * a.fps)

    def arm(q):
        T, _, _ = F.fk(q)
        return np.concatenate([(V @ T[L['frame']][:3, :3].T + T[L['frame']][:3, 3])[Fc]
                               for L, (V, Fc) in zip(links, meshes)])
    goal_tris = arm(qg)
    goal_col = np.c_[shade(goal_tris, GOAL), np.full(len(goal_tris), .15)]

    import imageio.v2 as imageio
    n = len(spec)
    fig = plt.figure(figsize=(4.2 * n, 4.6), dpi=110)
    axes = [fig.add_subplot(1, n, i + 1, projection='3d') for i in range(n)]
    fig.subplots_adjust(left=0, right=1, bottom=0, top=.86, wspace=0)
    writer = imageio.get_writer(a.out, fps=a.fps, codec='libx264', quality=8, macro_block_size=8)
    for fi in range(n_frames):
        t = fi / a.fps
        azim = -35 + 360 * fi / n_frames
        for ax, (f, label) in zip(axes, spec):
            ax.cla()
            q = traj[f][min(int(round(t / F.DT)), len(traj[f]) - 1)]
            A = arm(q)
            tris, cols = [tube, A, goal_tris], [tube_col, np.c_[shade(A, ROBOT), np.ones(len(A))], goal_col]
            ax.add_collection3d(Poly3DCollection(np.concatenate(tris), facecolors=np.concatenate(cols), edgecolor='none'))
            if f in spheres:
                ax.add_collection3d(Poly3DCollection(spheres[f], facecolor=TSPH + (.07,), edgecolor='none'))
            m = recs[f][k]
            done = t >= m['reached_s']
            ax.set_title(f"{label}\narrives {m['reached_s']:.2f} s, min gap "
                         + (f">= 20 mm" if m['min_gap_bound_mm'] >= 20 else f"{m['min_gap_bound_mm']:.1f} mm")
                         + ("\nREACHED" if done else f"\nt = {t:.2f} s"), fontsize=9, y=.98,
                         color=('green' if done else 'black'))
            ax.set_xlim(-.40, .50); ax.set_ylim(-.40, .40); ax.set_zlim(.15, .95)
            ax.set_box_aspect((1., .9, .9)); ax.view_init(18, azim); ax.set_axis_off()
        fig.canvas.draw()
        writer.append_data(np.asarray(fig.canvas.buffer_rgba())[..., :3])
        if fi % 20 == 0:
            print(f'frame {fi}/{n_frames}', flush=True)
    writer.close()
    print(a.out)


if __name__ == '__main__':
    main()
