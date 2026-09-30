"""Replay the heaviest steps of saved 6-mm sampled-data trajectories and solve each QP with the dense
least-distance solver and with Clarabel. python prototype_3d/qp_solver_benchmark.py 28 3"""
import os, sys, json, time, glob
for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'): os.environ[k] = '1'
os.environ['SUMMED_REFINE_DEPTH'] = '0'
sys.path.insert(0, '/home/user/t-ro-2027/tro/prototype_3d')
import numpy as np
import franka3d as F, summed as S, summed3d as S3, sampled_data as SD, certificate_upgrades as CU
R = '/home/user/t-ro-2027/tro/results'
links, obst, info = F.build(cache_path=R + '/fine_native_6mm_trial/cache_franka_6mm.npz')
S3.prep_all(links, obst.values()); SD.prepare(links, obst); CU.configure('sampled')
out = []
for trial in [int(a) for a in sys.argv[1:]]:
    d = np.load(f'{R}/certificate_upgrades/sampled_zero_6mm/franka_{trial:02d}.npz')
    Q, rows = d['q'], d['rows']
    kmax = int(np.argmax(rows)); ks = range(max(1, kmax - 6), min(len(rows), kmax + 4))
    zd = zs = None
    for k in ks:
        q, u_prev = Q[k], (Q[k] - Q[k - 1]) / F.DT
        st = SD.State(); st.u = u_prev
        qg = np.asarray(json.load(open('/home/user/t-ro-2027/tro/prototype_3d/franka_trials.json'))[trial][1])
        u_nom = np.clip(F.KQ * (qg - q), -F.QD_MAX, F.QD_MAX)
        t0 = time.perf_counter(); A, T, C, nb, hl, E, Wc = SD.franka_rows(q, links, obst, st); t_rows = time.perf_counter() - t0
        if Wc is not None and Wc.shape[1] == 0: Wc = None
        S.QP_SOLVER = 'nnls'; t0 = time.perf_counter(); ud, sd_, zd = S.solve(u_nom, A, T, C, E, F.BOX_G, F.BOX_H, zd, Wc); t_d = time.perf_counter() - t0
        S.QP_SOLVER = 'clarabel'; st0 = dict(S.QP_STATS); t0 = time.perf_counter(); us, ss, zs = S.solve(u_nom, A, T, C, E, F.BOX_G, F.BOX_H, zs, Wc); t_s = time.perf_counter() - t0
        path = 'sparse' if S.QP_STATS['sparse'] > st0['sparse'] else 'dense-fallback'
        rec = dict(trial=trial, step=k, rows=len(C), boxes=nb, t_rows_ms=1e3 * t_rows, t_dense_ms=1e3 * t_d, t_sparse_ms=1e3 * t_s,
                   path=path, du=float(np.abs(ud - us).max()), slack_dense=float(sd_))
        out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open(R + '/certificate_upgrades/qp_solver_benchmark.json', 'w'), indent=1)
