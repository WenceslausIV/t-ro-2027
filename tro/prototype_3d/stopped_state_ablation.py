"""At the stopped states of ours (12 mm, rows at the samples, results/certificate_upgrades/free), ablate the QP:
remove the velocity remainder (T = 0), or raise the constant terms by 5 cm, and measure the speed toward the
nominal direction. Run from tro/: python prototype_3d/stopped_state_ablation.py"""
import os, sys, json, glob
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['SUMMED_REFINE_DEPTH'] = '0'
sys.path.insert(0, 'prototype_3d')
import numpy as np, scipy.sparse as sp
import certificate_upgrades as CU, franka3d as F, summed as S, summed3d as S3
CU.configure('free'); S.QP_SOLVER = 'daqp'
links, obst, _ = F.build(); S3.prep_all(links, obst.values())
trials = json.load(open('prototype_3d/franka_trials.json'))
for f in sorted(glob.glob('results/certificate_upgrades/free/franka_*.npz')):
    i = int(f[-6:-4]); q = np.load(f)['q']
    if json.load(open(f[:-4] + '.json'))['metrics']['reached_s'] is not None or np.linalg.norm(q[-1] - q[-101]) > .02:
        continue
    qk, qg = q[-1], np.asarray(trials[i][1])
    un = np.clip(F.KQ * (qg - qk), -F.QD_MAX, F.QD_MAX)
    A, T, C, nb, hl, E, Wc = S3.franka_rows(qk, links, obst)
    prog = lambda u: float(u @ un / np.linalg.norm(un))          # speed toward the nominal direction [rad/s]
    u1, s1, _ = S.solve(un, A, T, C, E, F.BOX_G, F.BOX_H, None, Wc)
    T0 = sp.csr_matrix(T.shape) if sp.issparse(T) else np.zeros_like(T)
    u2, s2, _ = S.solve(un, A, T0, C, E, F.BOX_G, F.BOX_H, None, Wc)
    Cp = np.maximum(C, 0.) + .05 * F.GAMMA                         # constant terms raised by 5 cm (value-bound ablation)
    u3, s3, _ = S.solve(un, A, T, Cp, E, F.BOX_G, F.BOX_H, None, Wc)
    ref = []
    for depth in (1, 2):                                          # fix 1: halve near-contact boxes (online refinement)
        S.REFINE_DEPTH = depth
        A_, T_, C_, nb_, hl_, E_, Wc_ = S3.franka_rows(qk, links, obst)
        u4, s4, _ = S.solve(un, A_, T_, C_, E_, F.BOX_G, F.BOX_H, None, Wc_)
        ref.append(f'refine x{depth} {prog(u4):.3f} ({len(C_)} rows, slack {s4:.1e})')
    S.REFINE_DEPTH = 0
    print(f'trial {i:2d} rows {len(C):5d} | progress full {prog(u1):.3f}  no velocity remainder {prog(u2):.3f}  '
          f'+5cm constant terms {prog(u3):.3f}  ' + '  '.join(ref) + f'  (nominal {np.linalg.norm(un):.2f})', flush=True)
