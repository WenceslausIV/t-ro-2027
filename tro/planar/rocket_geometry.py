from pathlib import Path
"""Flame-free rocket silhouette traced from the user's reference illustration.

The dense polygon is the physical ground truth; all decorative patches are
clipped to that same silhouette. Heading zero points the nose along +x.
"""
import numpy as np
from matplotlib.patches import Polygon, Circle
import cspace_sdf_cbf_compare as C


def cubic(start, segments, samples=16):
    pieces=[];a=np.array(start,dtype=float)
    for c1,c2,end in segments:
        b=np.array([a,c1,c2,end],dtype=float)
        pieces.append(C.bern3(np.arange(samples)/samples)@b)
        a=np.array(end,dtype=float)
    return np.vstack(pieces)


OUTLINE=cubic((101,101),[
    ((131,94),(165,95),(182,98)),
    ((239,101),(267,126),(299,172)),
    ((321,202),(339,229),(353,250)),
    ((379,251),(409,249),(428,260)),
    ((443,268),(456,302),(475,329)),
    ((480,338),(472,338),(464,332)),
    ((449,322),(428,326),(415,340)),
    ((397,359),(382,380),(366,401)),
    ((350,402),(336,413),(329,429)),
    ((323,444),(329,458),(336,469)),
    ((342,481),(330,478),(321,471)),
    ((301,458),(272,440),(260,425)),
    ((251,414),(250,378),(249,355)),
    ((205,327),(164,302),(135,270)),
    ((110,242),(99,213),(96,183)),
    ((94,155),(98,124),(101,101)),
])
BODY=cubic((101,101),[
    ((131,94),(165,95),(182,98)),
    ((239,101),(267,126),(299,172)),
    ((332,217),(374,279),(415,340)),
    ((396,361),(381,382),(366,401)),
    ((302,388),(187,325),(135,270)),
    ((110,242),(99,213),(96,183)),
    ((94,155),(98,124),(101,101)),
])
AXES=np.array([[-1.,-1.],[1.,-1.]])/np.sqrt(2.)
RAW=OUTLINE@AXES.T
CENTER=(RAW.min(axis=0)+RAW.max(axis=0))/2
RADIUS=.30    # circumradius of the physical rocket [m]
SCALE=RADIUS/np.linalg.norm(RAW-CENTER,axis=1).max()


def local(points):
    return (np.asarray(points)@AXES.T-CENTER)*SCALE


def setup():
    from cspace_cbf_5robots import fit_curve,certify_curve_encloses
    gt=local(OUTLINE)
    area=np.sum(gt[:,0]*np.roll(gt[:,1],-1)-gt[:,1]*np.roll(gt[:,0],-1))
    if area<0:
        gt=gt[::-1]
    # Uniform arclength samples make the boundary fit insensitive to tracing density.
    closed=np.vstack([gt,gt[:1]])
    arc=np.r_[0,np.linalg.norm(np.diff(closed,axis=0),axis=1).cumsum()]
    q=np.linspace(0,arc[-1],1024,endpoint=False)
    points=np.column_stack([np.interp(q,arc,closed[:,j]) for j in range(2)])
    for offset in (.002,.003,.004,.006,.008,.010):
        ctrl=fit_curve(points,64,offset,n_fit=512)
        _,ok=certify_curve_encloses(ctrl,gt,max_depth=9)
        print(f'Rocket enclosure: offset {offset*1000:.1f} mm, passed {ok}',flush=True)
        if ok:
            robot=C.SmoothShape(ctrl,n_dense=768)
            return robot,gt,dict(name='reference rocket without exhaust',controls=64,
                fitting_offset_m=offset,enclosure_ok=True,physical_radius_m=float(np.linalg.norm(gt,axis=1).max()),
                fitted_radius_m=float(robot.rho),length_m=float(np.ptp(gt[:,0])),width_m=float(np.ptp(gt[:,1])))
    raise RuntimeError('Rocket body enclosure failed')


def draw_plain(ax,pose,gt,robot,alpha=1.,zorder=4,face='0.6',edge='#1f77b4',lw=.7):
    """Geometry only, paper convention: ground truth filled gray, certified fitted boundary colored."""
    ax.add_patch(Polygon(gt@C.rot(pose[2]).T+pose[:2],facecolor=face,edgecolor='none',alpha=alpha,zorder=zorder))
    boundary=robot.world(pose)
    boundary=np.vstack([boundary,boundary[:1]])
    ax.plot(boundary[:,0],boundary[:,1],color=edge,lw=lw,alpha=alpha,zorder=zorder+.1)


def draw(ax,pose,gt,robot,alpha=1.,zorder=4):
    def world(p):
        return local(p)@C.rot(pose[2]).T+pose[:2]
    physical=gt@C.rot(pose[2]).T+pose[:2]
    silhouette=Polygon(physical,facecolor='#f04c32',edgecolor='none',alpha=alpha,zorder=zorder)
    ax.add_patch(silhouette)
    parts=[(BODY,'#f5fbfb'),
           (np.array([(101,101),(182,98),(96,183)]),'#ff5936'),
           (np.array([(307,305),(407,401),(361,370)]),'#ff5936')]
    for poly,color in parts:
        patch=Polygon(world(poly),facecolor=color,edgecolor='none',alpha=alpha,zorder=zorder+.1)
        patch.set_clip_path(silhouette);ax.add_patch(patch)
    center=world(np.array([[203.,205.]]))[0]
    for radius,color in ((64,'#f04c32'),(52,'#07566b'),(46,'#04b2c8')):
        patch=Circle(center,radius*SCALE,facecolor=color,edgecolor='none',alpha=alpha,zorder=zorder+.2)
        patch.set_clip_path(silhouette);ax.add_patch(patch)
    boundary=robot.world(pose)
    boundary=np.vstack([boundary,boundary[:1]])
    ax.plot(boundary[:,0],boundary[:,1],color='#17425a',lw=.65,alpha=alpha,zorder=zorder+.3)


if __name__=='__main__':
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    robot,gt,meta=setup()
    print(meta)
    fig,ax=plt.subplots(figsize=(4,4))
    draw(ax,np.array([0.,0.,3*np.pi/4]),gt,robot)
    ax.set(xlim=(-.32,.32),ylim=(-.32,.32),aspect='equal');ax.axis('off')
    fig.savefig(str(Path(__file__).resolve().parents[1] / 'results' / 'rocket_geometry_preview.png'),dpi=200,bbox_inches='tight')
