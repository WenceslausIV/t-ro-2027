"""Evaluate the native-patch default, preserving every historical experiment.

python prototype_3d/native_patch_evaluation.py
"""
import os
for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['SUMMED_REFINE_DEPTH'] = '0'
import json
import argparse
from pathlib import Path
import numpy as np
import torch
import certificate_ablation as CA
import dock_cover as D
import franka3d as F
import summed as S
import summed3d as S3

OUT = CA.ROOT / 'results/native_patch_evaluation'


def metadata(scene, fields):
    return dict(scene=scene, cover_policy='native SDF patches, s=h', refinement_depth=0,
                row_mode='vertex', multiplier='one', dt_s=.01,
                geometry_sha256=CA.digest(fields), threads=1,
                timing_note='Background activity not excluded; exploratory timings.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial-range', nargs=2, type=int, default=[0,30], metavar=('START','STOP'))
    args=parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(1)
    S.REFINE_DEPTH, S.ROW_MODE, S.MULT_MODE = 0, 'vertex', 'one'
    if not (OUT/'dock.json').exists():
        gs,gc = D.station_body(),D.craft_body()
        fs,fc = D.body_field(gs,.0125,.15),D.body_field(gc,.01,.06)
        ms,mc = F.hessian_majorant(fs),F.hessian_majorant(fc)
        ls,lc = D.tight_level_2d(fs,ms,gs),D.tight_level_2d(fc,mc,gc)
        pc,side = D.cover_2d(fc,lc)
        assert side == fc.h
        fn = D.make_summed_fn(fs,ms,ls,S.prepare_body(fc,lc,pc,side,2))
        shapes = [type('Shape',(),dict(rho=float(np.linalg.norm(g,axis=1).max())))() for g in (gs,gc)]
        print('native docking:',len(pc),'boxes',flush=True)
        log = D.simulate_team([gs,gc],shapes,{(0,1):fn},D.DOCK2_START,D.DOCK2_GOAL,'summed',
                             steps=900,v_max=np.array([0.,1.]),w_max=np.array([0.,2.]),record=True)
        meta = metadata('dock',[fs.W,fc.W,[ls,lc],gs,gc])
        meta.update(side_mm=1000*side,cover_boxes=len(pc),levels_mm=[1000*ls,1000*lc])
        CA.save(OUT/'dock.json',log,meta,planar=True)
        p = log['final'][1]
        tip = p[:2]+D.NOSE_TIP_X*np.array([np.cos(p[2]),np.sin(p[2])])
        record=json.loads((OUT/'dock.json').read_text())
        record['metrics'].update(tip_depth_m=float(D.MOUTH_X-tip[0]),
                                 seats=bool(D.MOUTH_X-tip[0]>=D.SEAT_DEPTH),
                                 final_gap_mm=1000*D.C.true_distance(gs,gc@D.C.rot(p[2]).T+p[:2]))
        (OUT/'dock.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    links,obstacles,info=F.build()
    S3.prep_all(links,obstacles.values())
    assert all(abs(L['side']-L['sb']['f'].h)<1e-12 for L in links)
    reach=max(L['r'] for L in links)
    assert all(B['_M3_reach']>=reach for B in obstacles.values())
    arrays=[]
    for L in links: arrays += [L['sb']['f'].W,[L['l']],L['V'],L['F']]
    for B in obstacles.values():arrays += [B['f'].W,[B['l']]]
    meta=metadata('Franka 30 fixed trials',arrays)
    meta['setup']=info
    trials=json.loads((Path(F.HERE)/'franka_trials.json').read_text())[:30]
    print('native Franka covers:',[(L['frame'],len(L['PC']),L['side']) for L in links],flush=True)
    for i,(q0,qg,_) in enumerate(trials):
        if not args.trial_range[0] <= i < args.trial_range[1]:continue
        path=OUT/f'franka_{i:02d}.json'
        if path.exists():continue
        print('running native Franka',i,flush=True)
        log=F.simulate(q0,np.asarray(qg),links,obstacles,'ours',steps=1000,record=True)
        # Include the final state in the distance audit as well.
        final_gap=F.real_gap(log['q'][-1],links,obstacles)
        CA.save(path,log,dict(meta,trial_zero_based=i))
        record=json.loads(path.read_text())
        record['metrics']['final_gap_bound_mm']=1000*final_gap if np.isfinite(final_gap) else None
        path.write_text(json.dumps(record,indent=2),encoding='utf-8')
    if not all((OUT/f'franka_{i:02d}.json').exists() for i in range(len(trials))):
        print('Requested trial range complete; aggregate awaits other ranges.',flush=True)
        return
    import native_patch_report
    native_patch_report.main()


if __name__=='__main__':main()
