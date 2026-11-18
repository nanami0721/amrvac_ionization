#!/usr/bin/env python3
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
root,out=map(Path,sys.argv[1:3]);out.mkdir(parents=True,exist_ok=True)
rows=json.loads((root/'results.json').read_text());final=[r for r in rows if r['snapshot']==2]
by={(r['case'],r['placement'],r['direction'],r['scale']):r for r in final}
fig,axes=plt.subplots(1,3,figsize=(12.8,3.9),layout='constrained')
colors={1:'#137e91',-1:'#b45338'}
for direction in [1,-1]:
 for case in [0,1]:
  rr=[by[(case,1,direction,s)] for s in [1,2,3]]
  axes[case].plot([16,32,48],[r['roi_L1'] for r in rr],marker='o' if direction==1 else 's',ls='-' if direction==1 else '--',color=colors[direction],label='Forward (out of Yin)' if direction==1 else 'Reverse (into Yin)')
  control=by[(case,0,direction,2)];axes[case].scatter([32],[control['roi_L1']],marker='x',s=70,color=colors[direction],zorder=6)
 rr=[by[(1,1,direction,s)] for s in [1,2,3]]
 axes[2].plot([16,32,48],[r['div_defect_rms'] for r in rr],marker='o' if direction==1 else 's',ls='-' if direction==1 else '--',color=colors[direction])
for ax in axes:
 ax.set_xscale('log');ax.set_yscale('log');ax.set_xticks([16,32,48],['16','32','48']);ax.set_xlabel('Radial cells per patch');ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(alpha=.2,which='both');ax.spines[['top','right']].set_visible(False)
axes[0].set_title('Entropy pulse');axes[1].set_title('Force-free magnetic pulse');axes[2].set_title('Magnetic divergence defect')
axes[0].set_ylabel('Local relative L1 error');axes[1].set_ylabel('Local relative L1 error');axes[2].set_ylabel(r'RMS of $D_h(B_{num})-D_h(B_{exact})$')
axes[0].legend(fontsize=8);axes[1].text(.03,.04,'×: interior-region control',transform=axes[1].transAxes,fontsize=9)
fig.suptitle('Yin–Yang interface transport, t = 0.6 • identical mesh family and scheme')
fig.savefig(out/'convergence.png',dpi=180);fig.savefig(out/'convergence.pdf');plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(9,3.8),layout='constrained')
for i,metric in enumerate(['div_rms','seam_div_rms']):
 for direction in [1,-1]:
  rr=[by[(1,1,direction,s)] for s in [1,2,3]]
  axes[i].plot([16,32,48],[r[metric] for r in rr],marker='o' if direction==1 else 's',ls='-' if direction==1 else '--',color=colors[direction],label='Forward' if direction==1 else 'Reverse')
 if i==0:
  rr=[by[(1,1,1,s)] for s in [1,2,3]]
  axes[i].plot([16,32,48],[r['div_exact_sample_rms'] for r in rr],':',marker='x',color='#777777',label='Sampled analytic field (forward)')
 axes[i].set(xlabel='Radial cells per patch',ylabel='Diagnostic divergence RMS',title='Whole angular domain' if i==0 else 'Both sides of seam bands',xscale='log',yscale='log');axes[i].set_xticks([16,32,48],['16','32','48']);axes[i].xaxis.set_minor_formatter(NullFormatter());axes[i].grid(alpha=.2,which='both');axes[i].spines[['top','right']].set_visible(False)
axes[0].legend(fontsize=8);fig.suptitle('Raw divergence retained; angular boundary cells included')
fig.savefig(out/'divergence.png',dpi=180);fig.savefig(out/'divergence.pdf');plt.close(fig)
