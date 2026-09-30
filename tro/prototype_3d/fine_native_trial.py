"""Refit 6-mm link SDFs, use native cells as covers, and time one fixed Franka trial.

python prototype_3d/fine_native_trial.py
The first historically successful trial that stalled with 12-mm native patches is trial 2.
Obstacle fields are unchanged. No cover subdivision or online refinement is used.
"""
import os
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
    os.environ[key]='1'
os.environ['SUMMED_REFINE_DEPTH']='0'
import json
import platform
import time
from pathlib import Path
import numpy as np
import franka3d as F
import summed as S
import summed3d as S3

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/fine_native_6mm_trial'
TRIAL=2
H=.006


def stats(seconds):
    a=1000*np.asarray(seconds)
    return dict(median_ms=float(np.median(a)),p95_ms=float(np.percentile(a,95)),
                p99_ms=float(np.percentile(a,99)),max_ms=float(a.max()),mean_ms=float(a.mean()),
                over_10ms_steps=int(np.count_nonzero(a>10)),
                over_10ms_percent=float(100*np.mean(a>10)))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    S.REFINE_DEPTH,S.ROW_MODE,S.MULT_MODE=0,'vertex','one'
    old=np.load(F.CACHE,allow_pickle=True)['d'].item()
    fresh=[]
    fit_times=[]
    for i in F.LINK_IDS:
        p=OUT/f'link_{i}_6mm.npz'
        if p.exists():
            saved=np.load(p,allow_pickle=True)
            L=saved['link'].item(); elapsed=float(saved['setup_s'])
        else:
            print(f'Fitting and certifying link {i}: native SDF cell = 6 mm',flush=True)
            t=time.perf_counter()
            L=F.build_link(i,h=H)
            elapsed=time.perf_counter()-t
            np.savez_compressed(p,link=np.array(L,dtype=object),setup_s=elapsed)
        assert abs(L['f']['h']-H)<1e-12 and abs(L['side']-H)<1e-12
        fresh.append(L); fit_times.append(elapsed)
        print(f'link {i}: {len(L["PC"])} native patches, level {1000*L["l"]:.3f} mm, offline {elapsed:.2f} s',flush=True)
    cache=OUT/'cache_franka_6mm.npz'
    np.savez_compressed(cache,d=np.array(dict(links=fresh,obst=old['obst'],t_setup=sum(fit_times)),dtype=object))
    links,obst,info=F.build(cache_path=cache)
    t=time.perf_counter()
    S3.prep_all(links,obst.values())
    prep_s=time.perf_counter()-t
    assert all(abs(L['side']-L['sb']['f'].h)<1e-12 for L in links)
    trials=json.loads((Path(F.HERE)/'franka_trials.json').read_text())
    q0,qg,_=trials[TRIAL]
    q0,qg=np.asarray(q0),np.asarray(qg)
    print(f'Trial {TRIAL}: 6-mm refitted SDFs / native covers, no refinement; 10-ms budget',flush=True)
    # Prime numerical code before measuring the controller. No trajectory states are skipped.
    un=np.clip(F.KQ*(qg-q0),-F.QD_MAX,F.QD_MAX)
    Ar,Tr,Cr,nb,hl,Er=S3.franka_rows(q0,links,obst)
    S.solve(un,Ar,Tr,Cr,Er,F.BOX_G,F.BOX_H,None)
    q=q0.copy(); z=None; reached=None
    poses=[q.copy()]; us=[]; uns=[]; ts=[]; trs=[]; tqs=[]; rows=[]; boxes=[]; hs=[]; slacks=[]
    for k in range(1000):
        t0=time.perf_counter()
        un=np.clip(F.KQ*(qg-q),-F.QD_MAX,F.QD_MAX)
        Ar,Tr,Cr,nb,hl,Er=S3.franka_rows(q,links,obst)
        t1=time.perf_counter()
        u,slack,z=S.solve(un,Ar,Tr,Cr,Er,F.BOX_G,F.BOX_H,z)
        t2=time.perf_counter()
        ts.append(t2-t0); trs.append(t1-t0); tqs.append(t2-t1)
        us.append(u.copy());uns.append(un);rows.append(len(Cr));boxes.append(nb);hs.append(hl);slacks.append(slack)
        q=q+F.DT*u;poses.append(q.copy())
        if k%100==0:
            print(f'step {k}: rows {len(Cr)}, h lower bound {hl*1000:.3f} mm, goal error {np.linalg.norm(qg-q):.3f}',flush=True)
        if np.linalg.norm(qg-q)<.05:
            reached=(k+1)*F.DT
            break
    # Save timing and motion before the independent geometry audit.
    np.savez_compressed(OUT/'trajectory.npz',q=poses,u=us,u_nom=uns,time_s=ts,
                        assembly_s=trs,qp_s=tqs,rows=rows,boxes=boxes,h=hs,slack=slacks)
    print('Controller complete:',json.dumps(dict(reached_s=reached,filter=stats(ts))),flush=True)
    assert max(L['gt_r'].max() for L in links)<.005
    gaps=[]
    for k,q in enumerate(poses):
        gaps.append(min(.045,F.real_gap(q,links,obst)))
        if k%100==0:print(f'mesh audit {k}/{len(poses)-1}: running minimum {1000*min(gaps):.3f} mm',flush=True)
    np.save(OUT/'mesh_gap_lower_bounds_m.npy',gaps)
    historical=json.loads((Path(F.HERE)/'franka_results.json').read_text())['trials'][TRIAL]
    coarse=json.loads((ROOT/f'results/native_patch_evaluation/franka_{TRIAL:02d}.json').read_text())['metrics']
    meta=dict(trial_zero_based=TRIAL,selection='First historical success that stalled with the 12-mm native cover; selected before running.',
              link_sdf_cell_mm=6,training_spacing_mm=3,cover_side_mm=6,online_refinement_depth=0,
              obstacle_sdf_cell_mm=20,obstacle_fields_unchanged=True,refit=True,knot_insertion=False,
              row_mode='vertex',multiplier='one',dt_s=F.DT,threads=1,python=platform.python_version(),
              platform=platform.platform(),offline_link_fit_and_certify_s=sum(fit_times),offline_barrier_prep_s=prep_s,
              setup=info,timing_scope='Nominal command, FK, pruning, field evaluation, constraint assembly, QP. Excludes offline setup, mesh audit, rendering. One evaluation process; unrelated background activity not excluded.')
    result=dict(metadata=meta,steps=len(ts),reached_s=reached,final_goal_error=float(np.linalg.norm(qg-poses[-1])),
                cover_boxes=sum(len(L['PC']) for L in links),max_rows=max(rows),max_active_boxes=max(boxes),
                slack_steps=int(np.count_nonzero(np.asarray(slacks)>0)),min_h_mm=1000*min(hs),
                min_saved_mesh_gap_bound_mm=1000*min(gaps),saved_states_with_nonpositive_bound=int(np.count_nonzero(np.asarray(gaps)<=0)),
                filter=stats(ts),assembly=stats(trs),qp=stats(tqs),
                active_filter=stats(np.asarray(ts)[np.asarray(rows)>0]) if max(rows)>0 else None,
                historical_subdivided=dict(reached_s=historical['reached'],min_gap_mm=historical['min_gap_mm'],slack_steps=historical['slack_steps']),
                historical_12mm_native=coarse)
    (OUT/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    a=result['filter'];active=result['active_filter']
    report=f'''# Refit 6-mm native SDF patches: one Franka trial

Trial {TRIAL} was selected before execution: the first trial that succeeded with the historical
subdivided controller and stalled with 12-mm native patches. All eight link SDFs were **refitted**
with 6-mm cells and 3-mm training spacing using the existing distance-target routine; levels
were recertified on the meshes. Cover boxes are the native 6-mm cells. No knot insertion,
additional cover subdivision or online refinement was used. Obstacle SDFs remain at 20 mm.

| Setting | Goal arrival | Slack steps |
|---|---:|---:|
| Historical 12-mm SDF / 12-mm native cover | {coarse['reached_s']} | {coarse['slack_steps']} |
| Historical 12-mm SDF / 6-mm cover + optional refinement | {historical['reached']} s | {historical['slack_steps']} |
| Refitted 6-mm SDF / 6-mm native cover | {reached} s | {result['slack_steps']} |

Native cover: {result['cover_boxes']} boxes; maximum {result['max_active_boxes']} active boxes,
{result['max_rows']} rows. Minimum saved-state mesh-distance lower bound: {1000*min(gaps):.3f} mm.
The audit includes the final state, covers mesh triangles, caps geometrically pruned distances
at 45 mm, and does not certify between-step motion. Double precision, no outward rounding.

| Controller timing | ms |
|---|---:|
| Median | {a['median_ms']:.3f} |
| p95 | {a['p95_ms']:.3f} |
| p99 | {a['p99_ms']:.3f} |
| Maximum | {a['max_ms']:.3f} |
| Mean | {a['mean_ms']:.3f} |

10-ms budget exceeded on {a['over_10ms_steps']}/{len(ts)} steps ({a['over_10ms_percent']:.2f}%).
Active-step median/p95: {active['median_ms']:.3f}/{active['p95_ms']:.3f} ms.
This is a measured deadline check in Python, not a hard-real-time guarantee.
One CPU thread, one evaluation process, other background activity not excluded.
Filter timing includes nominal command, FK, pruning, spline evaluations, constraint assembly
and QP; collision audit and rendering run separately. Offline link fitting/certification:
{sum(fit_times):.2f} s; additional barrier preparation: {prep_s:.2f} s.

Reproduce: `python prototype_3d/fine_native_trial.py` (reuses this directory's fine-field cache).
The historical cache and manuscript are unchanged. One selected trial does not estimate a success rate.
'''
    (OUT/'README.md').write_text(report,encoding='utf-8')
    print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':main()
