"""Fixed-field, fixed-cover two-robot comparison; outputs only results/two_robot_cover_comparison.

python prototype_3d/two_robot_cover_comparison.py
python prototype_3d/two_robot_cover_comparison.py --render-only
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
import torch
import cspace_experiments as E
from cspace_cbf_5robots import random_shape
from dock_cover import body_field, cover_2d, tight_level_2d, make_summed_fn
import franka3d as F
import summed as SM

OUT = ROOT / 'results' / 'two_robot_cover_comparison'
DT, HORIZON = .01, 12.
STARTS = np.array([[-1.5, 0., .2], [.18, -1.45, 1.1]])
GOALS = np.array([[1.5, .15, -.7], [-.15, 1.5, 2.4]])
COLORS = ['#D47728', '#177FB5']
SIDES = {'patch': .04, 'fine': .005}
TITLES = {'patch': 'Patch as cover: 40 mm', 'fine': 'Subdivided cover: 5 mm'}
CORNERS = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])


def simulate(GT, fn, name):
    """Same controller and QP as simulate_team, run for a common fixed horizon."""
    x = STARTS.copy()
    G, h = E.box_rows(2, 'si', 1., 2.)
    radii = np.array([np.linalg.norm(g, axis=1).max() for g in GT])
    z = None
    data = {key: [] for key in ('traj', 'u', 'u_nom', 'time_ms', 'rows', 'boxes', 'h_low', 'slack')}
    arrival = None
    for k in range(round(HORIZON / DT)):
        tic = time.perf_counter()
        d = GOALS[:, :2] - x[:, :2]
        v = 1.5 * d
        v /= np.maximum(1., np.linalg.norm(v, axis=1))[:, None]
        nominal = np.c_[v, np.clip(2. * E.wrap(GOALS[:, 2] - x[:, 2]), -2., 2.)].ravel()
        result = fn(x[0], x[1]) if fn is not None else None
        if result is None:
            A, T, C, twist = np.zeros((0, 6)), np.zeros((0, 0)), np.zeros(0), np.zeros((0, 6))
            boxes, low = 0, np.nan
        else:
            A, T, C, boxes, low, twist = result
        u, slack, z = SM.solve(nominal, A, T, C, twist, G, h, z)
        elapsed = (time.perf_counter() - tic) * 1000
        values = (x.copy(), u.copy(), nominal, elapsed, len(C), boxes, low, slack)
        for key, value in zip(data, values):
            data[key].append(value)
        x += DT * u.reshape(2, 3)
        if arrival is None and np.all(np.linalg.norm(x[:, :2] - GOALS[:, :2], axis=1) < .05) and np.all(np.abs(E.wrap(x[:, 2] - GOALS[:, 2])) < .05):
            arrival = (k + 1) * DT
        if k % 200 == 0:
            print(f'{name}: t={k*DT:.1f}, rows={len(C)}, lower={low:.6f}, slack={slack:.3g}', flush=True)
    data = {k: np.asarray(v) for k, v in data.items()}
    data['final'] = x
    # Geometry audit is outside the controller timer and includes the final state.
    states = np.concatenate([data['traj'], x[None]])
    def polygon_gap(pose):
        world = [g @ E.rot(p[2]).T + p[:2] for g, p in zip(GT, pose)]
        return E.C.true_distance(*world)
    gaps, collision_states = [], 0
    for pose in states:
        signed_gap = polygon_gap(pose)
        gaps.append(max(0.,signed_gap))
        collision_states += int(signed_gap <= 0)
    # Hausdorff motion bound certifies each linear pose-interpolation interval
    # when max(endpoint gaps) exceeds the total possible relative displacement.
    # Otherwise recursively split; exact polygon distance and intersection at nodes.
    unresolved = 0
    interval_min = np.inf
    def audit_interval(a, b, da, db, depth=0):
        nonlocal unresolved
        motion = sum(np.linalg.norm(b[i, :2]-a[i, :2]) + radii[i]*abs(b[i, 2]-a[i, 2]) for i in range(2))
        lower = max(da, db) - motion
        if lower > 0:
            return lower
        if da == 0 or db == 0:
            return 0.
        if depth == 15:
            unresolved += 1
            return 0.
        mid = (a+b)/2
        dm = polygon_gap(mid)
        if dm <= 0:
            return 0.
        return min(audit_interval(a, mid, da, dm, depth+1), audit_interval(mid, b, dm, db, depth+1))
    for k in range(len(states)-1):
        interval_min = min(interval_min, audit_interval(states[k], states[k+1], gaps[k], gaps[k+1]))
    data['gap_mm'] = np.array(gaps) * 1000
    np.savez_compressed(OUT / f'{name}.npz', **data)
    active = data['rows'] > 0
    timing = data['time_ms']
    stats = dict(goal_time_s=arrival, collision_states=collision_states,
                 min_saved_gap_mm=min(gaps)*1000, interpolated_gap_lower_bound_mm=interval_min*1000,
                 unresolved_intervals=unresolved,
                 slack_steps=int(np.count_nonzero(data['slack'] > 0)),
                 max_slack=float(data['slack'].max()), max_rows=int(data['rows'].max()),
                 max_active_boxes=int(data['boxes'].max()), active_steps=int(active.sum()),
                 filter_median_ms=float(np.median(timing)), filter_p95_ms=float(np.percentile(timing,95)),
                 active_filter_median_ms=float(np.median(timing[active])) if active.any() else None,
                 active_filter_p95_ms=float(np.percentile(timing[active],95)) if active.any() else None,
                 input_modification=float(DT*np.linalg.norm(data['u']-data['u_nom'],axis=1).sum()),
                 path_lengths_m=np.linalg.norm(np.diff(states[:, :, :2],axis=0),axis=2).sum(axis=0).tolist(),
                 final_position_errors_m=np.linalg.norm(x[:,:2]-GOALS[:,:2],axis=1).tolist())
    print(name, json.dumps(stats), flush=True)
    return stats


def contours(f, level):
    axes = [np.arange(f.lo[i], f.lo[i]+f.K[i]*f.h, .002) for i in range(2)]
    X,Y = np.meshgrid(*axes)
    v = f.eval(np.c_[X.ravel(),Y.ravel(),np.zeros(X.size)]).reshape(X.shape)
    fig,ax = plt.subplots()
    paths = ax.contour(X,Y,v,levels=[level]).allsegs[0]
    plt.close(fig)
    return max(paths,key=len)


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(1)
    SM.REFINE_DEPTH = 0
    SM.ROW_MODE, SM.MULT_MODE = 'vertex', 'one'
    rng = np.random.default_rng(3)
    all_gt = [random_shape(rng) for _ in range(5)]
    GT = [all_gt[0],all_gt[3]]
    fields, majors, levels, curves = [], [], [], []
    for g in GT:
        f = body_field(g,.04,.1,training_spacing=.005)
        m = F.hessian_majorant(f)
        level = tight_level_2d(f,m,g)
        fields.append(f); majors.append(m); levels.append(level); curves.append(contours(f,level))
    np.savez_compressed(OUT/'geometry.npz',g0=GT[0],g1=GT[1],c0=curves[0],c1=curves[1])
    result = dict(dt_s=DT,horizon_s=HORIZON,patch_side_m=.04,training_spacing_m=.005,
                  levels_mm=(np.array(levels)*1000).tolist(), patch_counts=[int(np.prod(f.K[:2])) for f in fields],
                  starts=STARTS.tolist(),goals=GOALS.tolist(),refinement_depth=0,
                  row_mode='vertex',multiplier='one',activation_m=.05,
                  timing_note='One CPU thread; background activity not excluded. Fixed 12-s horizon, filter only; excludes geometry audit and drawing.',
                  audit_note='Exact polygon geometry: vertex-to-edge distances plus containment and edge-crossing tests at all states; motion bounds with recursive bisection for linear pose interpolation, double precision (no outward rounding).',
                  trials={})
    for name,side in SIDES.items():
        centers,actual_side = cover_2d(fields[1],levels[1],side)
        assert abs(actual_side-side)<1e-12
        body = SM.prepare_body(fields[1],levels[1],centers,actual_side,2)
        # make_summed_fn caches a radius-dependent majorant on the target field.
        # Reset it for each cover radius; both runs retain identical spline weights.
        if hasattr(fields[0],'_M3'):
            delattr(fields[0],'_M3')
        fn = make_summed_fn(fields[0],majors[0],levels[0],body)
        np.savez_compressed(OUT/f'{name}_cover.npz',centers=centers,side=side)
        print(f'{name}: {len(centers)} boxes, side={side}',flush=True)
        result['trials'][name] = dict(cover_boxes=len(centers),cover_side_mm=side*1000,**simulate(GT,fn,name))
        (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    result['trials']['unfiltered'] = simulate(GT,None,'unfiltered')
    (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')


def render():
    report = json.loads((OUT/'results.json').read_text())
    geometry = np.load(OUT/'geometry.npz')
    logs = {name:np.load(OUT/f'{name}.npz') for name in SIDES}
    covers = {name:np.load(OUT/f'{name}_cover.npz') for name in SIDES}
    fig,axs = plt.subplots(1,2,figsize=(12,6.3),dpi=100)
    fig.subplots_adjust(left=.04,right=.98,bottom=.12,top=.85,wspace=.12)
    fig.suptitle('Two moving robots | identical 40 mm SDF patches and controller',fontsize=15,y=.98)
    fig.text(.5,.925,'Only cover size changes. Automatic refinement OFF. Orange: robot A; blue: covered robot B.',ha='center',fontsize=10)
    fig.text(.5,.025,'Same simulation clock | dt = 10 ms | Dashed outlines: goals | Squares: entire surface cover of robot B',ha='center',fontsize=10)
    artists = {}
    for ax,(name,side) in zip(axs,SIDES.items()):
        ax.set_title(f'{TITLES[name]}  |  {report["trials"][name]["cover_boxes"]} boxes',fontsize=12)
        ax.set(xlim=(-2.05,2.05),ylim=(-2.05,2.05),aspect='equal',xlabel='x [m]',ylabel='y [m]')
        ax.grid(alpha=.15)
        paths,polys,lines = [],[],[]
        for i,color in enumerate(COLORS):
            g = geometry[f'g{i}']
            ax.add_patch(Polygon(g@E.rot(GOALS[i,2]).T+GOALS[i,:2],fill=False,edgecolor=color,ls='--',lw=.9,alpha=.5))
            ax.text(*GOALS[i,:2],f'G{i+1}',ha='center',va='center',color=color,fontsize=9)
            path, = ax.plot([],[],color=color,alpha=.5,lw=1.4)
            poly = Polygon(g,facecolor=color,edgecolor=color,alpha=.22,lw=.5,zorder=3)
            ax.add_patch(poly)
            line, = ax.plot([],[],color=color,lw=1,zorder=4)
            paths.append(path); polys.append(poly); lines.append(line)
        boxes = PolyCollection([],facecolors='none',edgecolors=COLORS[1],linewidths=.45 if name=='patch' else .25,alpha=.7,zorder=2)
        ax.add_collection(boxes)
        text = ax.text(.02,.98,'',transform=ax.transAxes,va='top',fontsize=9,family='monospace',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'),zorder=6)
        artists[name] = paths,polys,lines,boxes,text
    def update(k):
        for name in SIDES:
            log,cov = logs[name],covers[name]
            pose = log['traj'][k]
            paths,polys,lines,boxes,text = artists[name]
            for i in range(2):
                p = pose[i]
                world = geometry[f'g{i}']@E.rot(p[2]).T+p[:2]
                curve = geometry[f'c{i}']@E.rot(p[2]).T+p[:2]
                polys[i].set_xy(world); lines[i].set_data(curve[:,0],curve[:,1])
                paths[i].set_data(log['traj'][:k+1,i,0],log['traj'][:k+1,i,1])
            local = cov['centers'][:,None,:]+CORNERS[None]*float(cov['side'])
            boxes.set_verts(local@E.rot(pose[1,2]).T+pose[1,:2])
            goal = report['trials'][name]['goal_time_s']
            status = 'GOALS REACHED' if goal is not None and k*DT>=goal else 'MOVING'
            if status == 'MOVING' and np.linalg.norm(log['u'][k]) < .01:
                status = 'SLOW / STOPPED'
            if k==len(log['traj'])-1 and goal is None:
                status='GOALS NOT REACHED'
            text.set_text(f't = {k*DT:5.2f} s   {status}\ngap = {log["gap_mm"][k]:6.1f} mm   rows = {log["rows"][k]:4d}')
        return []
    frames = list(range(0,len(logs['patch']['traj']),8))+[len(logs['patch']['traj'])-1]
    anim = FuncAnimation(fig,update,frames=frames,interval=80,blit=False)
    anim.save(OUT/'comparison.gif',writer=PillowWriter(fps=12.5))
    contact = int(np.argmin(logs['fine']['gap_mm'][:-1]))
    update(contact)
    fig.savefig(OUT/'comparison.png',dpi=170)
    for ax,name in zip(axs,SIDES):
        pose = logs[name]['traj'][contact]
        world = [geometry[f'g{i}']@E.rot(pose[i,2]).T+pose[i,:2] for i in range(2)]
        _,pa,pb,_ = E.closest_pair(*(torch.as_tensor(p,device=E.DEV) for p in world))
        center = (pa+pb)/2
        ax.set_xlim(center[0]-.24,center[0]+.24)
        ax.set_ylim(center[1]-.24,center[1]+.24)
    fig.suptitle('Closest-approach detail | same 0.48 m window size in both panels',fontsize=15,y=.98)
    fig.savefig(OUT/'contact_zoom.png',dpi=200)
    plt.close(fig)
    fig,axs = plt.subplots(3,1,figsize=(9,7),sharex=True)
    for name,color in [('patch','#8F4DA8'),('fine','#158C83')]:
        log=logs[name]; ts=np.arange(len(log['time_ms']))*DT
        axs[0].plot(np.arange(len(log['gap_mm']))*DT,log['gap_mm'],color=color,label=TITLES[name])
        axs[1].plot(ts,log['rows'],color=color)
        axs[2].plot(ts,log['time_ms'],color=color,alpha=.65,lw=.6)
    axs[0].set(ylabel='Polygon gap [mm]',ylim=(0,150));axs[0].legend()
    axs[1].set(ylabel='QP barrier rows')
    axs[2].set(ylabel='Filter time [ms]',xlabel='Simulation time [s]',yscale='log')
    for ax in axs:ax.grid(alpha=.2)
    fig.suptitle('Same 12 s horizon; one CPU thread; background activity not excluded')
    fig.tight_layout();fig.savefig(OUT/'performance.png',dpi=160);plt.close(fig)
    rows=[]
    for name in SIDES:
        s=report['trials'][name]
        rows.append(f'| {TITLES[name]} | {s["goal_time_s"]} | {s["min_saved_gap_mm"]:.3f} | {s["interpolated_gap_lower_bound_mm"]:.3f} | {s["max_rows"]} | {s["filter_median_ms"]:.3f}/{s["filter_p95_ms"]:.3f} | {s["active_filter_median_ms"]:.3f}/{s["active_filter_p95_ms"]:.3f} | {s["input_modification"]:.3f} | {s["slack_steps"]} |')
    readme='''# Two-robot fixed-cover comparison

Same ground-truth shapes (seed 3, shape indices 0 and 3), 40 mm spline patches,
5 mm training spacing, certified levels, starts/goals, input bounds, gamma=5,
50 mm activation threshold, unit lifted vertex certificate, dt=10 ms, 12 s horizon.
Only the cover size changes: 40 mm native patches vs 5 mm subdivisions.
Automatic refinement is OFF. One direction per pair: robot B's surface against A's field.
The GIF shows all B cover boxes; only active boxes contribute QP rows.

| Cover | Both goals [s] | Min saved gap [mm] | Interpolated lower bound [mm] | Max rows | Filter median/p95 [ms] | Active filter median/p95 [ms] | Input modification J | Slack steps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+f'''

Unfiltered baseline: {report['trials']['unfiltered']['collision_states']} of 1201 saved states intersect.
Goal tolerance: position < 5 cm AND orientation error < 0.05 rad, both robots.
J = integral ||u-u_nom||_2 dt; mixes translational and angular coordinates identically in both runs.
Timings cover nominal command, constraint assembly, and QP, excluding collision audit and drawing.
Single CPU thread; unrelated background activity was not excluded. Active-step timings use each run's own active steps.
This is one predetermined crossing scene, not a success-rate benchmark or proof of universal improvement.
Polygon vertex-to-edge distances, containment, and edge-crossing tests check saved states. Interpolated lower bounds use
rigid-motion bounds and recursive bisection for the simulated linear pose interpolation in double precision;
they are not a sampled-data controller theorem or an outward-rounded numerical certificate.

Reproduce: `python prototype_3d/two_robot_cover_comparison.py`
Render saved results: `python prototype_3d/two_robot_cover_comparison.py --render-only`
'''
    (OUT/'README.md').write_text(readme,encoding='utf-8')
    print('Saved GIF, PNG figures, README and raw logs to',OUT,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--render-only',action='store_true')
    args=parser.parse_args()
    if not args.render_only:run()
    render()
