"""Refit the Franka link SDFs with a given cell size; the native cells are the cover patches (box = patch).

    python prototype_3d/refit_native_fields.py 18 24
Writes results/native_patch_sizes/<h>mm/link_<i>.npz and cache_franka_<h>mm.npz (obstacle fields unchanged).
Uses franka3d.build_link, the routine of the 12-mm and 6-mm fields: distance targets at spacing h/2 and the
tightest certified level on the link mesh. Box side = native cell = h; no subdivision, no merging.
"""
import os
import sys
import time
from pathlib import Path

for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np

import franka3d as F

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results' / 'native_patch_sizes'


def cache_path(mm):
    return OUT / f'{mm}mm' / f'cache_franka_{mm}mm.npz'


def refit(mm):
    h = mm / 1000
    d = OUT / f'{mm}mm'
    d.mkdir(parents=True, exist_ok=True)
    old = np.load(F.CACHE, allow_pickle=True)['d'].item()
    links, times = [], []
    for i in F.LINK_IDS:
        p = d / f'link_{i}.npz'
        if p.exists():
            s = np.load(p, allow_pickle=True)
            L, dt = s['link'].item(), float(s['setup_s'])
        else:
            t = time.perf_counter()
            L = F.build_link(i, h=h)
            dt = time.perf_counter() - t
            np.savez_compressed(p, link=np.array(L, dtype=object), setup_s=dt)
        assert abs(L['f']['h'] - h) < 1e-12 and abs(L['side'] - h) < 1e-12
        links.append(L); times.append(dt)
        print(f'{mm} mm link {i}: {len(L["PC"])} native patches, level {1e3 * L["l"]:.3f} mm, {dt:.1f} s', flush=True)
    np.savez_compressed(cache_path(mm), d=np.array(dict(links=links, obst=old['obst'], t_setup=sum(times)),
                                                   dtype=object))


if __name__ == '__main__':
    for mm in sys.argv[1:]:
        refit(int(mm))
