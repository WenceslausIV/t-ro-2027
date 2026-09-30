"""Curved maze (our own designed passage): vector geometry, one SE(2) field, and CBF escape.

python planar/reference_maze_escape.py [--layout-only] [--rebuild]
"""
import argparse
import hashlib
import heapq
import json
import time
from pathlib import Path
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import cspace_sdf_cbf_compare as C
import cspace_experiments as E
from maze_escape import simulate
import maze_reference_geometry as G
import rocket_geometry as Rocket

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results'/'reference_maze_rocket'
FIG = ROOT/'tro'/'figs'
DOMAIN, KXY, KTH = 6., 160, 48
RASTER, N_TRAIN_TH, N_TRAIN_XY = .02, 96, 300    # training targets of the field
PLANNER_MARGIN = .02
LOOKAHEAD = 1.2            # pure-pursuit look-ahead along the route [m]
WAYPOINT_CUT = .10         # waypoint segments may enter the walls by this much [m] (negative: stay this far from them)


def plan(robot,walls,half_width):
    """Offline A* supplies the nominal route; it is not a property of the CBF."""
    res=.03
    n=2*int(round(DOMAIN/res))+1
    g=(np.arange(n)-n//2)*res
    occupied=np.zeros((n,n),dtype=bool)
    for wall in walls:
        occupied |= C._raster_polygon(wall.dense,n,res).astype(bool)
    distance=ndimage.distance_transform_edt(~occupied)*res
    X,Y=np.meshgrid(g,g,indexing='ij')
    bounds=G.world([(34,410),(328,66)])
    free=(distance > half_width+PLANNER_MARGIN) & (X>bounds[0,0]) & (X<bounds[1,0]) & (Y>bounds[0,1]) & (Y<bounds[1,1])
    start,end=G.world([(65,83),(280,402)])
    index=lambda p: tuple(np.round(p/res+n//2).astype(int))
    src,dst=index(start),index(end)
    assert free[src] and free[dst], 'Start or exit is not in the inflated free space'
    cost={src:0.}; parent={}; queue=[(0.,src)]
    while queue:
        _,u=heapq.heappop(queue)
        if u==dst:
            break
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)):
            v=(u[0]+dx,u[1]+dy)
            if not (0<=v[0]<n and 0<=v[1]<n and free[v]):
                continue
            if dx and dy and not (free[u[0]+dx,u[1]] and free[u[0],u[1]+dy]):
                continue
            weight=np.hypot(dx,dy)*res*(1+.08/max(distance[v]-half_width,.025))
            new=cost[u]+weight
            if new<cost.get(v,np.inf):
                cost[v]=new; parent[v]=u
                heuristic=np.hypot(v[0]-dst[0],v[1]-dst[1])*res
                heapq.heappush(queue,(new+heuristic,v))
    if dst not in parent:
        raise RuntimeError('No route through the traced maze at this robot size')
    route=[dst]
    while route[-1]!=src:
        route.append(parent[route[-1]])
    return (np.array(route[::-1])-n//2)*res


def dead_end_fillers(walls,path,res=.02,keep_factor=1.25,keep_extra=.04):
    """Close the dead ends: free space farther from the route than the local corridor half-width
    (x keep_factor, + keep_extra) becomes smooth filler walls that overlap the original walls."""
    from scipy.spatial import cKDTree
    from skimage import measure
    n=2*int(round(DOMAIN/res))+1
    g=(np.arange(n)-n//2)*res
    occupied=np.zeros((n,n),dtype=bool)
    for wall in walls:
        occupied |= C._raster_polygon(wall.dense,n,res).astype(bool)
    distance=ndimage.distance_transform_edt(~occupied)*res
    X,Y=np.meshgrid(g,g,indexing='ij')
    bounds=G.world([(34,410),(328,66)])
    inside=(X>bounds[0,0]) & (X<bounds[1,0]) & (Y>bounds[0,1]) & (Y<bounds[1,1])
    free=~occupied & inside
    pix=np.argwhere(free)
    d,idx=cKDTree(path).query((pix-n//2)*res)
    pidx=np.clip(np.round(path/res+n//2).astype(int),0,n-1)
    r=distance[pidx[:,0],pidx[:,1]]
    keep=np.zeros((n,n),dtype=bool)
    keep[tuple(pix[d<=keep_factor*r[idx]+keep_extra].T)]=True
    pocket=ndimage.binary_opening(free & ~keep,iterations=2)
    lab,count=ndimage.label(pocket)
    sizes=ndimage.sum(pocket,lab,range(1,count+1))*res*res
    fill=np.isin(lab,1+np.flatnonzero(sizes>.05))
    fill=ndimage.binary_dilation(fill,iterations=5) & (occupied | fill)      # overlap into the walls only
    smooth=ndimage.gaussian_filter(np.pad(fill.astype(float),8),4.)
    fillers=[]
    for c in measure.find_contours(smooth,.5):
        world=(c-8-n//2)*res
        closed=np.vstack([world,world[:1]])
        area=.5*abs(np.sum(closed[:-1,0]*closed[1:,1]-closed[1:,0]*closed[:-1,1]))
        if area>.05:
            fillers.append(G.ContourWall(world))
    return fillers


def smooth_corridor_walls(walls,path,res=.02,sigma_center=.3,sigma_width=.6,width_scale=1.0,width_min=0.):
    """Walls = maze frame minus ONE smooth tube: its centerline is the route smoothed along
    arclength, its half-width the local corridor half-width of the traced maze smoothed along
    arclength. Dead ends disappear and every wall boundary is a smooth offset curve."""
    from scipy.spatial import cKDTree
    from skimage import measure
    n=2*int(round(DOMAIN/res))+1
    g=(np.arange(n)-n//2)*res
    occupied=np.zeros((n,n),dtype=bool)
    for wall in walls:
        occupied |= C._raster_polygon(wall.dense,n,res).astype(bool)
    distance=ndimage.distance_transform_edt(~occupied)*res
    # resample the route uniformly and smooth it; extend it straight out of the frame at both ends
    arc=np.r_[0,np.linalg.norm(np.diff(path,axis=0),axis=1).cumsum()]
    ds=.02; q=np.arange(0,arc[-1],ds)
    xy=np.column_stack([np.interp(q,arc,path[:,i]) for i in range(2)])
    pidx=np.clip(np.round(xy/res+n//2).astype(int),0,n-1)
    width=np.maximum(width_scale*ndimage.gaussian_filter1d(distance[pidx[:,0],pidx[:,1]],sigma_width/ds,mode='nearest'),width_min)
    center=ndimage.gaussian_filter1d(xy,sigma_center/ds,axis=0,mode='nearest')
    head=np.column_stack([np.full(40,center[0,0]),center[0,1]+ds*np.arange(40,0,-1)])
    tail=np.column_stack([np.full(40,center[-1,0]),center[-1,1]-ds*np.arange(1,41)])
    center=np.vstack([head,center,tail]); width=np.r_[np.full(40,width[0]),width,np.full(40,width[-1])]
    X,Y=np.meshgrid(g,g,indexing='ij')
    pts=np.column_stack([X.ravel(),Y.ravel()])
    d,idx=cKDTree(center).query(pts)
    corridor=(d<=width[idx]).reshape(n,n)
    lo,hi=G.world([(33,411),(329,65)])                       # outer frame of the maze
    frame=(X>=lo[0]) & (X<=hi[0]) & (Y>=lo[1]) & (Y<=hi[1])
    soft=ndimage.gaussian_filter(corridor.astype(float),5.)
    wall_field=np.where(frame,1.-soft,0.)
    # round off thin wedges (e.g., where a wall meets the frame at an acute angle)
    rr=np.hypot(*np.mgrid[-6:7,-6:7])<=6
    wall_field=ndimage.gaussian_filter(ndimage.binary_opening(wall_field>.5,structure=rr).astype(float),5.)
    out=[]
    for c in measure.find_contours(np.pad(wall_field,2),.5):
        world=(c-2-n//2)*res
        closed=np.vstack([world,world[:1]])
        area=.5*abs(np.sum(closed[:-1,0]*closed[1:,1]-closed[1:,0]*closed[:-1,1]))
        if area>.05:
            out.append(G.ContourWall(world,spacing=.03))
    return out,center[40:-40]          # walls and the smooth centerline (inside the frame)


FRAME = np.array([[-4.44,-3.69],[4.44,5.19]])        # square 8.88 x 8.88 m
# Our own winding passage (not traced from any image): the entrance drops vertically at the top left, three
# horizontal passes (each with one smooth wiggle that vanishes to first order at its ends) are joined by
# semicircular hairpins of radius 1.1 m, and the exit leaves vertically at the bottom right. The centerline is
# tangent-continuous, its radius of curvature is at least ~0.9 m, larger than every half-width, so both wall
# boundaries are smooth offset curves, and the passage crosses the frame at right angles.
TURN_R, PASS_Y, PASS_X, WIGGLE = 1.1, (3.3, 1.1, -1.1), 2.2, .2
HALF_WIDTH = (.47,.04,4.1)          # half-width = a + b sin(2 pi s / c) along the arclength s [m]
CORNER_RADIUS = .06                 # smooth max radius at every wall corner (entrance, exit, frame)


def snake_centerline(ds=.004):
    def line(p,q):
        n=max(2,int(np.linalg.norm(np.subtract(q,p))/ds)); t=np.linspace(0,1,n)[:,None]
        return np.asarray(p)+t*(np.subtract(q,p))
    def arc(c,r,a0,a1):
        n=max(2,int(abs(a1-a0)*r/ds)); a=np.linspace(a0,a1,n)
        return np.column_stack([c[0]+r*np.cos(a),c[1]+r*np.sin(a)])
    def wiggly(x0,x1,y):
        n=int(abs(x1-x0)/ds); u=np.linspace(0,1,n)
        return np.column_stack([x0+(x1-x0)*u,y+WIGGLE*np.sin(2*np.pi*u)*np.sin(np.pi*u)])
    R,X=TURN_R,PASS_X
    y1,y2,y3=PASS_Y
    parts=[line((-X-R,6.3),(-X-R,y1+R)), arc((-X,y1+R),R,np.pi,1.5*np.pi),       # entrance, turn into pass 1
           wiggly(-X,X,y1), arc((X,y1-R),R,.5*np.pi,-.5*np.pi),                   # pass 1, hairpin
           wiggly(X,-X,y2), arc((-X,y2-R),R,.5*np.pi,1.5*np.pi),                  # pass 2, hairpin
           wiggly(-X,X,y3), arc((X,y3-R),R,.5*np.pi,0.),                          # pass 3, turn down
           line((X+R,y3-R),(X+R,-6.3))]                                           # exit
    return np.vstack(parts)


CENTER_POINTS = snake_centerline(2)[[0,-1]]      # entrance and exit (for the labels)


def designed_corridor_walls(res=.01):
    """Walls = frame minus ONE smooth tube around a designed centerline. Both sets are signed-distance
    fields, combined with a polynomial smooth max, so every corner (also at the entrance and exit, where the
    tube leaves the frame) is rounded with radius ~CORNER_RADIUS; walls are the zero contour on a 1-cm grid."""
    from scipy.spatial import cKDTree
    from skimage import measure
    c=snake_centerline()
    s=np.r_[0,np.linalg.norm(np.diff(c,axis=0),axis=1).cumsum()]
    a,b,period=HALF_WIDTH
    w=a+b*np.sin(2*np.pi*s/period)
    lo,hi=FRAME
    xs=np.arange(lo[0]-.3,hi[0]+.3+1e-9,res); ys=np.arange(lo[1]-.3,hi[1]+.3+1e-9,res)
    X,Y=np.meshgrid(xs,ys,indexing='ij'); P=np.column_stack([X.ravel(),Y.ravel()])
    d,i=cKDTree(c).query(P)
    tube=d-w[i]                                                           # <0 inside the passage
    q=np.abs(P-(lo+hi)/2)-(hi-lo)/2
    frame=np.linalg.norm(np.maximum(q,0),axis=1)+np.minimum(q.max(axis=1),0)  # <0 inside the frame
    k=CORNER_RADIUS
    aa,bb=frame,-tube                                                     # wall = frame AND NOT tube
    hmix=np.maximum(k-np.abs(aa-bb),0)/k
    F=(np.maximum(aa,bb)+hmix**2*k/4).reshape(X.shape)                    # smooth max (rounded corners)
    out=[]
    for cont in measure.find_contours(F,0.):
        world=np.column_stack([np.interp(cont[:,0],np.arange(len(xs)),xs),np.interp(cont[:,1],np.arange(len(ys)),ys)])
        closed=np.vstack([world,world[:1]])
        area=.5*abs(np.sum(closed[:-1,0]*closed[1:,1]-closed[1:,0]*closed[:-1,1]))
        if area>.05:
            out.append(G.ContourWall(world,spacing=.03))
    inside=(c[:,0]>lo[0])&(c[:,0]<hi[0])&(c[:,1]>lo[1])&(c[:,1]<hi[1])
    route=c[inside][::20]
    return out,route


def few_waypoints(path,walls,cut=.10,res=.02):
    """As few goal waypoints as possible: from each kept route point, jump to the farthest route
    point whose straight segment enters the walls by at most `cut` (signed distance >= -cut), so
    the CBF, not the waypoints, keeps the rocket off the walls. The last waypoint is the exit."""
    n=2*int(round(DOMAIN/res))+1
    occupied=np.zeros((n,n),dtype=bool)
    for wall in walls:
        occupied |= C._raster_polygon(wall.dense,n,res).astype(bool)
    sd=(ndimage.distance_transform_edt(~occupied)-ndimage.distance_transform_edt(occupied))*res
    def ok(p,q):
        t=np.linspace(0,1,max(int(np.linalg.norm(q-p)/(res/2)),2))[:,None]
        idx=np.round((p+t*(q-p))/res+n//2).astype(int)
        return sd[idx[:,0],idx[:,1]].min()>=-cut
    keep=[0]
    while keep[-1]<len(path)-1:
        i=keep[-1]; j=len(path)-1
        while j>i+1 and not ok(path[i],path[j]):
            j-=1
        keep.append(j)
    xy=path[keep]
    seg=np.diff(xy,axis=0)
    heading=np.r_[np.arctan2(seg[0,1],seg[0,0]),np.arctan2(seg[:,1],seg[:,0])]
    return np.column_stack([xy,heading])


def field(robot,walls,rebuild):
    config=dict(domain=DOMAIN,cells=[KXY,KTH],robot=robot.ctrl.tolist(),
        walls=[w.bez.tolist() for w in walls],raster=RASTER,n_xy=N_TRAIN_XY,n_theta=N_TRAIN_TH)
    signature=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    cache,mp=OUT/'field.npz',OUT/'field.json'
    f=C.Field3D('bspline',DOMAIN,KXY,KTH)
    if cache.exists() and mp.exists() and not rebuild:
        meta=json.loads(mp.read_text())
        if meta['signature']==signature:
            f.set_W(np.load(cache)['W'])
            if all(meta[k] for k in ('contact_ok','band_ok','regular_ok')):
                return f,meta
    start=time.perf_counter()
    res=RASTER; D=DOMAIN+.1; n=2*int(round(D/res))+1
    g=(np.arange(n)-n//2)*res
    obs=np.zeros((n,n),dtype=float)
    for wall in walls:
        obs=np.maximum(obs,C._raster_polygon(wall.dense,n,res))
    th=np.arange(N_TRAIN_TH)*C.TWO_PI/N_TRAIN_TH
    slices=np.empty((len(th),n,n),dtype=np.float32)
    print('Fitting the reference maze field...',flush=True)
    for k,angle in enumerate(th):
        rb=C._raster_polygon(-robot.dense@C.rot(angle).T,n,res)
        occ=fftconvolve(obs,rb,mode='same')>.5
        slices[k]=(ndimage.distance_transform_edt(~occ)-ndimage.distance_transform_edt(occ))*res
    xy=np.linspace(-DOMAIN,DOMAIN,N_TRAIN_XY)
    f.fit(C.sample_slices(g,slices.astype(float),xy,xy).transpose(1,2,0),xy,xy,th)
    fitting=time.perf_counter()-start
    print(f'Field fitted: {fitting:.1f} s; checking contact and regularity...',flush=True)
    cert=C.Certifier(f,robot,E.Union(walls))
    level,ok,regions,tc=cert.certify(n_theta0=8,verbose=True,chunk=100000)
    band,bmin=cert.certify_band(level)
    regular,depth,remaining=cert.certify_regular(level,max_depth=12)
    meta=dict(signature=signature,level=float(level),fit_s=fitting,contact_s=tc,
        contact_ok=bool(ok),band_ok=bool(band),band_min=float(bmin),regular_ok=bool(regular),
        regular_depth=depth,regular_open_boxes=remaining,regions=regions,cells=[KXY,KXY,KTH],
        coefficients=f.n_coef,wall_chord_error_m=max(w.chord_error for w in walls))
    np.savez_compressed(cache,W=f.W); mp.write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2),flush=True)
    if not (ok and band and regular):
        raise RuntimeError('Incomplete certification; no paper results exported')
    return f,meta


def draw(robot,gt,walls,arrays,preview=False,fitted_walls=None,out_dir=None,fig_path=None):
    out=OUT if out_dir is None else out_dir
    fig_path=FIG/'rocket_maze.png' if fig_path is None else fig_path
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(4.9,5.15))
    LW=1.6                      # ~0.9 pt at 0.8 column width (the paper includes it at 0.8 of the column)
    # colors as in Fig. 1: static environment neutral (light-gray fill, dark-gray certified outline), the
    # robot in the one accent hue (light-blue fill, blue certified outline, blue path), the goal green
    WALL_FILL,WALL_EDGE,ROCKET_FILL,ROCKET_EDGE='0.85','0.35','#c6dbef','#1f77b4'
    for wall in walls:
        ax.add_patch(Polygon(wall.dense,facecolor=WALL_FILL,edgecolor='none',zorder=3))
    if fitted_walls is not None:
        for wall in fitted_walls:
            curve=np.vstack([wall.dense,wall.dense[:1]])
            ax.plot(curve[:,0],curve[:,1],color=WALL_EDGE,lw=LW,zorder=3.2)   # on top of the fill, like the rocket outlines
    poses=arrays['poses']
    arc=np.r_[0,np.linalg.norm(np.diff(poses[:,:2],axis=0),axis=1).cumsum()]
    if preview:
        ax.plot(poses[:,0],poses[:,1],color='#2477c5',lw=1.25,ls='--',zorder=3)
    else:                       # path opacity proportional to the arclength travelled (as in the five-robot figure)
        from matplotlib.collections import LineCollection
        from matplotlib.colors import to_rgba
        q=np.arange(0,arc[-1]+.01,.01).clip(max=arc[-1])
        pts=np.column_stack([np.interp(q,arc,poses[:,k]) for k in range(2)])
        seg=np.stack([pts[:-1],pts[1:]],axis=1)
        rgba=np.tile(to_rgba('#2a78d6'),(len(seg),1)); rgba[:,3]=(q[:-1]+q[1:])/(2*arc[-1])
        ax.add_collection(LineCollection(seg,colors=rgba,linewidths=2.9,capstyle='butt',zorder=3.5))
    ids=[] if preview else list(np.searchsorted(arc,np.arange(.9,arc[-1]-.5,1.7)))
    # one rocket seen at several times: its opacity follows the path's (faint early, strong late), with a
    # floor so that the first pose stays visible; the final pose lies under the goal disc and is omitted
    for i in [0]+ids:
        a=1. if preview else .15+.85*arc[i]/arc[-1]
        Rocket.draw_plain(ax,poses[i],gt,robot,alpha=a,zorder=4,face=ROCKET_FILL,edge=ROCKET_EDGE,lw=LW)
    # goal (end of the path, in the exit): flat exit-sign green disc marked G (the figure is planar)
    ex,rg=poses[-1,:2],.24
    ax.add_patch(Circle(ex,rg,facecolor='#00a651',edgecolor='none',zorder=6))
    ax.text(*ex,'G',color='white',ha='center',va='center',fontsize=11,weight='bold',zorder=7)
    pad=FRAME[0,0]+4.8                                # same margin below the frame as left of it
    ax.set(xlim=(-4.8,4.8),ylim=(FRAME[0,1]-pad,FRAME[1,1]+.35),aspect='equal',xlabel='$x$ [m]',ylabel='$y$ [m]')
    ax.spines[['top','right']].set_visible(False)
    ax.set_xticks(np.arange(-4,5,2)); ax.set_yticks(np.arange(-2,5,2))
    if not preview:
        fig.legend(handles=[Patch(facecolor=WALL_FILL,edgecolor=WALL_EDGE,label='Maze walls'),
            Line2D([],[],color='#2a78d6',lw=2.,label='Rocket trajectory (faint: early)')],
            loc='lower center',ncol=2,fontsize=8,frameon=False)
    fig.tight_layout(pad=.5,rect=(0,.08 if not preview else 0,1,1))
    if preview:
        fig.savefig(out/'layout_preview.png',dpi=200,bbox_inches='tight')
    else:
        for dest in (fig_path,out/'maze_escape.png'):
            fig.savefig(dest,dpi=250,bbox_inches='tight')
        fig.savefig(out/'rocket_maze.pdf',bbox_inches='tight')     # the paper uses the PNG
    plt.close(fig)
    if not preview:
        fig,ax=plt.subplots(figsize=(5.8,5.2))
        for wall,fit in zip(walls,fitted_walls):
            ax.add_patch(Polygon(wall.dense,facecolor='#b9bdc3',edgecolor='none'))
            curve=np.vstack([fit.dense,fit.dense[:1]])
            ax.plot(curve[:,0],curve[:,1],color='#c53131',lw=2.)
        ax.set(xlim=(-4.47,-4.35),ylim=(5.12,5.23),aspect='equal',xlabel='$x$ [m]',ylabel='$y$ [m]')
        ax.legend(handles=[Patch(facecolor='#b9bdc3',label='Physical wall'),
            Line2D([],[],color='#252525',lw=2,label='Enclosing cubic B-spline')],fontsize=8)
        fig.tight_layout(); fig.savefig(out/'tight_wall_corner.png',dpi=220); plt.close(fig)
        fig,axes=plt.subplots(2,1,figsize=(7,4.6),sharex=True)
        t=np.arange(len(poses))*.01
        axes[0].plot(t,arrays['gaps'].min(axis=1)*1000); axes[0].set_ylabel('Gap lower bound [mm]')
        axes[1].plot(t,arrays['barriers']*1000); axes[1].axhline(0,color='k',ls='--',lw=.7)
        axes[1].set(ylabel='$h$ [mm]',xlabel='Time [s]'); fig.tight_layout()
        fig.savefig(out/'clearance.png',dpi=160); plt.close(fig)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--layout-only',action='store_true')
    ap.add_argument('--rebuild',action='store_true'); ap.add_argument('--redraw',action='store_true')
    args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    robot,gt,robot_meta=Rocket.setup()
    np.savez_compressed(OUT/'robot_geometry.npz',gt=gt,ctrl=robot.ctrl)
    walls,path=designed_corridor_walls()                  # our own passage; path = its centerline
    print(f'Designed corridor: {len(walls)} wall regions',flush=True)
    fitted_walls=[]; wall_fits=[]
    for wall in walls:
        fitted,fit_meta=G.fit_wall(wall)
        fitted_walls.append(fitted); wall_fits.append(fit_meta)
    np.savez_compressed(OUT/'wall_fits.npz',**{f'{key}_{i}':getattr(w,key)
        for i,w in enumerate(fitted_walls) for key in ('ctrl','knots','bez','span_lengths')})
    if args.redraw:                                   # figure only, from the saved trajectory
        draw(robot,gt,walls,dict(np.load(OUT/'trajectory.npz')),fitted_walls=fitted_walls); return
    waypoints=few_waypoints(path,walls,cut=WAYPOINT_CUT)
    print(f'Route: {len(waypoints)} waypoints, robot radius {robot.rho:.3f} m',flush=True)
    if args.layout_only:
        draw(robot,gt,walls,dict(poses=waypoints),preview=True,fitted_walls=fitted_walls); return
    f,meta=field(robot,fitted_walls,args.rebuild)
    arrays,result=simulate(f,meta['level'],gt,walls,waypoints=waypoints,domain=DOMAIN,
        gap_function=G.physical_gaps,output_dir=OUT,cruise=True,steps=18000,switch_radius=.40,dyn='si')
    files=['reference_maze_escape.py','maze_reference_geometry.py','maze_escape.py','rocket_geometry.py']
    result.update(field=meta,robot=robot_meta,planner_margin=PLANNER_MARGIN,wall_count=len(walls),wall_fits=wall_fits,
        source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
        navigation='Dead ends closed; fewest goal waypoints whose segments enter walls by <= WAYPOINT_CUT; CBF filters controls.',
        geometry='Traced cubic physical walls enclosed by separately fitted periodic cubic B-splines.',
        gap_definition='Polygon distance minus the uniform wall-curve chord deviation bound.',
        verification='Float64 certificates and interval rigid-motion audit; no outward rounding.')
    (OUT/'results.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(OUT/'trajectory.npz',**arrays)
    draw(robot,gt,walls,arrays,fitted_walls=fitted_walls)
    metrics=dict(mazeTime=f"{result['reached_s']:.1f}",mazeGap=f"{1000*result['min_physical_gap_m']:.1f}",
        mazeLower=f"{1000*result['verified_interval_gap_lower_m']:.1f}",mazeLevel=f"{1000*meta['level']:.1f}",
        mazeMedian=f"{result['controller_median_ms']:.2f}",mazeInterventions=str(result['barrier_intervention_steps']))
    (ROOT/'tro'/'maze_results.tex').write_text('% Generated by reference_maze_escape.py\n'+
        ''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key,value in metrics.items()))


if __name__=='__main__':
    main()
