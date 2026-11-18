#!/usr/bin/env python3
"""Plot field lines traced through actual native Yin–Yang dipole snapshots.
Usage: python plot_dipole.py SOURCE RUN_DIR OUTPUT_DIR
Requires NumPy, SciPy, Matplotlib; run outside the source tree.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LogNorm
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import check_mhd
source,run,out=map(Path,sys.argv[1:4]);out.mkdir(parents=True,exist_ok=True)
util=check_mhd.reader(source)
tilt=np.deg2rad(20);m=np.array([np.sin(tilt),0,np.cos(tilt)])
rotation=np.array([[-1,0,0],[0,0,1],[0,1,0]])
interpolators=[];bounds=[];metrics={}

def exact(pos):
    r=np.linalg.norm(pos,axis=-1);n=pos/r[...,None]
    return .5*(3*np.sum(n*m,axis=-1)[...,None]*n-m)/r[...,None]**3

for patch,name in enumerate(['yin','yang']):
    h,w=check_mhd.load(util,run/f'yy_{name}0000.dat')
    dx=(h['xmax']-h['xmin'])/h['domain_nx']
    axes=[a+(np.arange(n)+.5)*d for a,n,d in zip(h['xmin'],h['domain_nx'],dx)]
    r,t,p=np.meshgrid(*axes,indexing='ij');st,ct,sp,cp=np.sin(t),np.cos(t),np.sin(p),np.cos(p)
    basis=np.stack([np.stack([st*cp,st*sp,ct],-1),np.stack([ct*cp,ct*sp,-st],-1),np.stack([-sp,cp,np.zeros_like(p)],-1)],-1)
    if patch:basis=np.einsum('ij,...jk->...ik',rotation,basis)
    pos=r[...,None]*basis[..., :,0]
    field=np.einsum('...ij,...j->...i',basis,w[...,5:8])
    metrics[name+'_initial_B_max_abs_error']=float(np.max(abs(field-exact(pos))))
    interpolators.append(RegularGridInterpolator(axes,field,bounds_error=False,fill_value=None))
    bounds.append([(a[0],a[-1]) for a in axes])
    ht,wt=check_mhd.load(util,run/f'yy_{name}0001.dat')
    pressure=(ht['gamma']-1)*(wt[...,4]-.5*np.sum(wt[...,1:4]**2,-1)/wt[...,0]-.5*np.sum(wt[...,5:8]**2,-1))
    metrics[name+'_five_step_time']=ht['time']
    metrics[name+'_five_step_B_relative_l2_change']=float(np.linalg.norm(wt[...,5:8]-w[...,5:8])/np.linalg.norm(w[...,5:8]))
    metrics[name+'_five_step_pressure_min']=float(pressure.min())

switches=np.zeros(2,dtype=int)
def sample(points):
    points=np.atleast_2d(points);coords=[];scores=[]
    for patch in range(2):
        local=points if patch==0 else points@rotation.T
        rr=np.linalg.norm(local,axis=1);tt=np.arctan2(np.linalg.norm(local[:,:2],axis=1),local[:,2]);pp=np.arctan2(local[:,1],local[:,0])
        coords.append(np.column_stack([rr,tt,pp]))
        scores.append(np.minimum.reduce([tt-bounds[patch][1][0],bounds[patch][1][1]-tt,pp-bounds[patch][2][0],bounds[patch][2][1]-pp]))
    choose=np.argmax(scores,axis=0)
    if np.any(np.max(scores,axis=0)<-1e-12):raise RuntimeError('point has no angular interior donor')
    value=np.empty_like(points)
    for patch in range(2):
        mask=choose==patch
        if np.any(mask):value[mask]=interpolators[patch](coords[patch][mask]);switches[patch]+=mask.sum()
    return value

# Seed in the magnetic frame. Field lines use interpolated snapshot vectors.
e1=np.array([np.cos(tilt),0,-np.sin(tilt)]);e2=np.array([0,1,0.])
seeds=[];directions=[]
for deg in [23,32,43,53,63,73]:
    for phi in np.deg2rad(np.arange(0,360,45)):
        th=np.deg2rad(deg)
        seeds.append(1.04*(np.cos(th)*m+np.sin(th)*(np.cos(phi)*e1+np.sin(phi)*e2)));directions.append(1)
        if deg<40:
            seeds.append(1.04*(-np.cos(th)*m+np.sin(th)*(np.cos(phi)*e1+np.sin(phi)*e2)));directions.append(-1)
x=np.array(seeds);directions=np.array(directions);paths=[[p.copy()] for p in x];active=np.ones(len(x),bool);step=.018
for iteration in range(1300):
    ids=np.flatnonzero(active)
    if not len(ids):break
    xx=x[ids];sgn=directions[ids,None]
    def rhs(q):
        field=sample(q);return sgn*field/np.linalg.norm(field,axis=1)[:,None]
    k1=rhs(xx);k2=rhs(xx+.5*step*k1);k3=rhs(xx+.5*step*k2);k4=rhs(xx+step*k3)
    new=xx+step*(k1+2*k2+2*k3+k4)/6
    for i,q in zip(ids,new):paths[i].append(q.copy())
    x[ids]=new;rad=np.linalg.norm(new,axis=1);active[ids]=(rad>=1.035)&(rad<=3.)
paths=[np.array(p) for p in paths];metrics['traced_lines']=len(paths);metrics['unfinished_lines']=int(active.sum());metrics['sample_count_per_patch']=switches.tolist()
# Check the dipole field-line invariant r/sin(theta_m)^2 along numerical traces.
invariants=[]
for p in paths:
    rr=np.linalg.norm(p,axis=1);ct=p@m/rr;inv=rr/(1-ct*ct)
    invariants.append(float(np.max(abs(inv/inv[0]-1))))
metrics['fieldline_invariant_max_relative_drift']=max(invariants)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titleweight':'bold'})
fig=plt.figure(figsize=(13.4,6.6),facecolor='#101a29')
ax=fig.add_axes([.025,.14,.53,.76],facecolor='#101a29');ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-3.15,3.15);ax.set_ylim(-3.15,3.15)
# Orthographic projection with explicit sphere occultation (no 3-D painter artifacts).
az,el=np.deg2rad(38),np.deg2rad(17)
view=np.array([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)])
right=np.array([-np.sin(az),np.cos(az),0]);up=np.cross(view,right)
q=np.linspace(-1,1,500);X,Y=np.meshgrid(q,q);inside=X*X+Y*Y<=1;Z=np.sqrt(np.maximum(0,1-X*X-Y*Y))
normal=X[...,None]*right+Y[...,None]*up+Z[...,None]*view
Br=10*np.sum(normal*m,axis=-1);cmap=plt.get_cmap('RdBu_r');rgba=cmap(Normalize(-10,10)(Br));rgba[...,:3]*=(.52+.48*Z)[...,None];rgba[...,3]=inside
ax.imshow(rgba,extent=[-1,1,-1,1],origin='lower',zorder=3)
ax.add_patch(Circle((0,0),1,fill=False,edgecolor='#c3d4e4',lw=.8,zorder=4))
for path in paths:
    a,b,depth=path@right,path@up,path@view
    visible=~((a*a+b*b<1)&(depth<0));p2=np.column_stack([a,b]);segments=np.stack([p2[:-1],p2[1:]],1)
    ax.add_collection(LineCollection(segments[visible[:-1]&visible[1:]],colors='#ffd276',linewidths=.85,alpha=.78,zorder=5))
# Dipole-axis markers and outer tracing limit.
ax.add_patch(Circle((0,0),3,fill=False,ls=(0,(3,5)),edgecolor='#577087',lw=.7,zorder=0))
ax.text(-2.95,-2.93,r'$3\,R_\odot$ tracing limit',color='#b4c9d9',fontsize=10)
axis2=np.array([m@right,m@up]);ax.annotate('',xy=1.43*axis2,xytext=1.05*axis2,arrowprops=dict(arrowstyle='->',color='white',lw=1.5));ax.text(*(1.57*axis2),'N',color='white',ha='center',va='center',weight='bold')
ax.set_title('Tilted solar dipole',color='white',fontsize=19,pad=15)
cax=fig.add_axes([.11,.105,.32,.016]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-10,10),cmap=cmap),cax=cax,orientation='horizontal');cb.set_label(r'Photospheric $B_r$ [G]  (analytic boundary)',color='#dce7ef',fontsize=10);cb.ax.tick_params(colors='#dce7ef',labelsize=9);cb.outline.set_edgecolor('#708394')

ax2=fig.add_axes([.615,.22,.31,.63],facecolor='#101a29');axis=np.linspace(-3,3,360);xx,zz=np.meshgrid(axis,axis);rr=np.hypot(xx,zz);mask=(rr>=1.04)&(rr<=2.99);pts=np.column_stack([xx[mask],np.zeros(mask.sum()),zz[mask]])
field=sample(pts);strength=np.full_like(xx,np.nan);strength[mask]=10*np.linalg.norm(field,axis=1)
im=ax2.pcolormesh(xx,zz,strength,norm=LogNorm(.18,10),cmap='magma',shading='auto',rasterized=True)
ax2.add_patch(Circle((0,0),1,facecolor='#e9eff4',edgecolor='white',lw=1.2,zorder=4));ax2.text(0,0,'Sun',ha='center',va='center',color='#36465a',fontsize=14,zorder=6)
# Streamplot on the numerical meridional sample; interior/exterior points masked.
bx=np.full_like(xx,np.nan);bz=bx.copy();bx[mask]=field[:,0];bz[mask]=field[:,2]
ax2.streamplot(axis,axis,np.ma.masked_invalid(bx),np.ma.masked_invalid(bz),color='#f7e9ca',density=.85,linewidth=.6,arrowsize=.7)
ax2.plot([0,m[0]*1.6],[0,m[2]*1.6],ls='--',color='#708092',lw=1.2,zorder=5)
ax2.set_aspect('equal');ax2.set_xlim(-3,3);ax2.set_ylim(-3,3);ax2.set_xlabel(r'$x/R_\odot$',color='#dce7ef');ax2.set_ylabel(r'$z/R_\odot$',color='#dce7ef');ax2.tick_params(colors='#dce7ef');ax2.set_title('Snapshot field in the x–z plane',color='white',fontsize=13,pad=12)
for spine in ax2.spines.values():spine.set_color('#71869a')
cax2=fig.add_axes([.94,.25,.012,.56]);cb2=fig.colorbar(im,cax=cax2);cb2.set_label(r'$|B|$ [G]',color='#dce7ef');cb2.ax.tick_params(colors='#dce7ef',labelsize=9);cb2.outline.set_edgecolor('#708394')
fig.text(.62,.055,'t = 0 • 32 × 32 × 96 cells per patch\nPotential field; no wind, gravity or heating',color='#c2d1dd',fontsize=10,linespacing=1.7)
fig.text(.03,.965,'AMRVAC / YIN–YANG DEMO',color='#8fa9bd',fontsize=10,weight='bold')
fig.savefig(out/'solar_dipole.png',dpi=190,facecolor=fig.get_facecolor());fig.savefig(out/'solar_dipole.pdf',facecolor=fig.get_facecolor());plt.close(fig)
(out/'validation.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(metrics,indent=2))
