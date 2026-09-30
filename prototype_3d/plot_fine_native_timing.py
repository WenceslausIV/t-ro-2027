"""Plot measured controller timing for the saved 6-mm native-patch trial."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent.parent/'results/fine_native_6mm_trial'
r=json.loads((OUT/'result.json').read_text())
z=np.load(OUT/'trajectory.npz')
t=np.arange(len(z['time_s']))*.01
fig,ax=plt.subplots(3,1,figsize=(9,7),sharex=True,constrained_layout=True)
ax[0].plot(t,1000*z['time_s'],color='#245a9b',lw=.8,label='Complete controller')
ax[0].axhline(10,color='#c34738',ls='--',label='100 Hz budget (10 ms)')
ax[0].set_ylabel('Compute time [ms]');ax[0].legend(loc='upper right')
ax[1].plot(t,z['rows'],color='#245a9b',lw=1)
ax[1].set_ylabel('CBF rows')
g=np.load(OUT/'mesh_gap_lower_bounds_m.npy')*1000
ax[2].plot(np.arange(len(g))*.01,g,color='#227755',lw=1)
ax[2].axhline(0,color='#c34738',lw=.7)
ax[2].set_ylabel('Mesh gap bound [mm]');ax[2].set_xlabel('Simulated time [s]')
for a in ax:a.grid(alpha=.2);a.spines[['top','right']].set_visible(False)
goal=f"Goal reached at {r['reached_s']:.2f} s" if r['reached_s'] is not None else 'Goal not reached within 10 s'
fig.suptitle('Refitted 6-mm SDF = 6-mm native cover; no refinement\n'
             f"{goal} | >10 ms: {r['filter']['over_10ms_percent']:.1f}% of steps",fontsize=12)
fig.savefig(OUT/'timing.png',dpi=180)
print(OUT/'timing.png')
