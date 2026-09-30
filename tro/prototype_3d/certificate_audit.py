"""Audit saved Franka runs: does every executed input satisfy the certificate rows built at its state?

    python prototype_3d/certificate_audit.py <results subfolder> {unit|free|sampled} {12mm|6mm}
Rebuilds the rows at every saved state (native patches, no refinement) and checks the executed input
u_k = (q_{k+1} - q_k) / dt exactly: A u + T |E u| + C >= 0, with, for optimized multipliers, some w in [0, W_MAX]
per box; zero inputs are skipped (always admissible). Rows with A = 0 are the sampled-data caps and are counted
separately. Writes results/certificate_upgrades/audit/<folder>.json.
"""
import os, sys, glob, json
os.environ['OMP_NUM_THREADS'] = '1'; os.environ.setdefault('SUMMED_REFINE_DEPTH', '0')
import numpy as np, scipy.sparse as sp
import certificate_upgrades as CU, franka3d as F, summed as S, summed3d as S3, sampled_data as SD
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
folder, variant, field = sys.argv[1], sys.argv[2], sys.argv[3]
CU.configure(variant)
if field == '6mm':
    links, obst, info = F.build(cache_path=ROOT / 'results/fine_native_6mm_trial/cache_franka_6mm.npz')
else:
    links, obst, info = F.build()
S3.prep_all(links, obst.values())
if variant == 'sampled':
    SD.prepare(links, obst)

def check(u, A, T, C, E, Wc):
    Td = T.toarray() if sp.issparse(T) else T
    al = A @ u + Td @ np.abs(E @ u) + C
    tol = 1e-9 * np.maximum(1., np.abs(C))
    capmask = ~np.any(A != 0, axis=1)
    if Wc is None or Wc.shape[1] == 0:
        bad = al < -tol
        return int((bad & ~capmask).sum()), int((bad & capmask).sum()), float(al.min()) if len(al) else 0.
    W = sp.coo_matrix(Wc); inW = np.zeros(len(C), bool); inW[W.row] = True
    bad = (al < -tol) & ~inW
    nb = W.shape[1]; lo = np.zeros(nb); hi = np.full(nb, S.W_MAX)
    pos, neg = W.data > 0, W.data < 0
    np.maximum.at(lo, W.col[pos], -al[W.row[pos]] / W.data[pos])
    np.minimum.at(hi, W.col[neg], -al[W.row[neg]] / W.data[neg])
    boxbad = int((lo > hi + 1e-9).sum())
    return int((bad & ~capmask).sum()) + boxbad, int((bad & capmask).sum()), float(np.min(hi - lo)) if nb else 0.

trials = json.loads((ROOT / 'prototype_3d/franka_trials.json').read_text())
out = {}
for f in sorted(glob.glob(str(ROOT / 'results' / folder / 'franka_*.npz'))):
    i = int(f[-6:-4]); q = np.load(f)['q']; rec = dict(steps=len(q) - 1, barrier_bad_steps=0, cap_bad_steps=0,
                                                     input_bad_steps=0, first_bad=None)
    up = np.zeros(F.DOF)
    for k in range(len(q) - 1):
        u = (q[k + 1] - q[k]) / F.DT
        if variant == 'sampled':
            st = SD.State(); st.u = up
            A, T, C, nb, hl, E, Wc = SD.franka_rows(q[k], links, obst, st)
        else:
            r = S3.franka_rows(q[k], links, obst); A, T, C, nb, hl, E = r[:6]; Wc = r[6] if len(r) == 7 else None
        up = u
        if not np.any(u) or len(C) == 0:
            continue
        b, c, _ = check(u, A, T, C, E, Wc)
        rec['barrier_bad_steps'] += b > 0; rec['cap_bad_steps'] += c > 0
        rec['input_bad_steps'] += bool(np.any(F.BOX_G @ u - F.BOX_H < -1e-9))
        if b > 0 and rec['first_bad'] is None:
            rec['first_bad'] = (k, b)
    out[i] = rec
    print(folder, i, rec, flush=True)
tot = {k: sum(r[k] for r in out.values()) for k in ('steps', 'barrier_bad_steps', 'cap_bad_steps', 'input_bad_steps')}
print('TOTAL', folder, tot, 'trials with barrier violations', [i for i, r in out.items() if r['barrier_bad_steps']], flush=True)
dst = ROOT / 'results' / 'certificate_upgrades' / 'audit'
dst.mkdir(parents=True, exist_ok=True)
(dst / (folder.replace('/', '_') + '.json')).write_text(json.dumps(dict(folder=folder, variant=variant, field=field,
                                                                     per_trial=out, total=tot), indent=1))
