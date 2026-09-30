"""Compare the physical maze with a tight enclosing workspace-SDF level set.

This is a geometry illustration, separate from the SE(2) escape controller.
Run: python maze_levelset_figure.py
"""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.path import Path as PlotPath
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

import cspace_sdf_cbf_compare as C
from cspace_cbf_5robots import sub_matrix_np
import maze_reference_geometry as G

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results'/'maze_workspace_levelset'
FIG = ROOT/'tro'/'figs'
DOMAIN, CELLS, RASTER, TOL = 6., 600, .0025, .0001


def fit(walls):
    # Constant angular basis: exactly a bicubic workspace field, with no robot.
    f = C.Field3D('bspline', DOMAIN, CELLS, 1)
    n = 2*round(DOMAIN/RASTER)+1
    occupied = np.zeros((n,n),dtype=bool)
    for wall in walls:
        occupied |= C._raster_polygon(wall.dense,n,RASTER).astype(bool)
    sdf = (ndimage.distance_transform_edt(~occupied)-ndimage.distance_transform_edt(occupied))*RASTER
    xy = np.linspace(-DOMAIN,DOMAIN,1801)
    X,Y = np.meshgrid((xy+DOMAIN)/RASTER,(xy+DOMAIN)/RASTER,indexing='ij')
    target = ndimage.map_coordinates(sdf,[X,Y],order=1)
    f.fit(target[:,:,None],xy,xy,np.array([0.]))
    return f


def outside_boxes(centers, radii, walls):
    """True only if a box misses both continuous physical wall regions."""
    outside = np.ones(len(centers),dtype=bool)
    for wall in walls:
        coarse=wall.coarse
        lower=coarse.mid_tree.query(centers)[0]-coarse.max_half_edge-coarse.chord_error
        separated=(~coarse.path.contains_points(centers)) & (lower>radii)
        near=np.flatnonzero((lower<=radii)&outside)
        points=centers[near]
        upper = wall.vertex_tree.query(points)[0]
        ids = wall.mid_tree.query_ball_point(points,upper+wall.max_half_edge+1e-10)
        distance = np.empty(len(points))
        for j,(p,edges) in enumerate(zip(points,ids)):
            a=wall.edge_a[edges]; v=wall.edge_b[edges]-a
            t=np.clip(np.sum((p-a)*v,axis=1)/np.sum(v*v,axis=1),0,1)
            distance[j]=np.linalg.norm(p-a-t[:,None]*v,axis=1).min()
        separated[near]=(~wall.path.contains_points(points)) & (distance>radii[near]+wall.chord_error)
        outside &= separated
    return outside


def refined_walls(walls):
    """Tighten the auxiliary polygon tubes for the whole-region certificate."""
    refined=[]
    for wall in walls:
        pieces=[];errors=[]
        for b in wall.bez:
            m2=np.linalg.norm(6*np.diff(b,n=2,axis=0),axis=1).max()
            n=max(4,int(np.ceil(np.sqrt(m2/(8e-7)))),
                  int(np.ceil(np.linalg.norm(np.diff(b,axis=0),axis=1).sum()/.02)))
            pieces.append(C.bern3(np.arange(n)/n)@b);errors.append(m2/(8*n*n))
        dense=np.vstack(pieces);ends=np.roll(dense,-1,axis=0)
        refined.append(SimpleNamespace(coarse=wall,path=PlotPath(dense),
            edge_a=dense,edge_b=ends,chord_error=max(errors),
            mid_tree=cKDTree((dense+ends)/2),vertex_tree=cKDTree(dense),
            max_half_edge=np.linalg.norm(ends-dense,axis=1).max()/2))
    return refined


def certify_interior(f,level,walls):
    """Check phi < level on the whole physical wall union, not only its boundary."""
    coeff = f.C[:,:,0,:,:,0]
    ids = np.argwhere(coeff.max(axis=(2,3)) >= level)
    lo = -DOMAIN+ids*f.bx.h
    widths = np.full((len(ids),2),f.bx.h)
    polys = coeff[ids[:,0],ids[:,1]]
    halves=sub_matrix_np(np.array([0.,.5]),np.array([.5,1.]))
    checked=0
    for depth in range(18):
        checked += len(polys)
        if not len(polys):
            return dict(interior_ok=True,interior_depth=depth,interior_boxes=checked)
        keep=~outside_boxes(lo+widths/2,np.linalg.norm(widths,axis=1)/2,walls)
        polys,lo,widths=polys[keep],lo[keep],widths[keep]
        children=[]; positions=[]
        for ix in range(2):
            for iy in range(2):
                sub=np.einsum('ip,npq,jq->nij',halves[ix],polys,halves[iy])
                select=sub.max(axis=(1,2))>=level
                children.append(sub[select])
                positions.append((lo+widths*np.array([ix,iy])/2)[select])
        polys=np.concatenate(children)
        lo=np.concatenate(positions)
        widths=np.full((len(polys),2),f.bx.h/2**(depth+1))
        print(f'Interior check depth {depth}: {len(polys)} unresolved boxes',flush=True)
    raise RuntimeError('Whole-wall enclosure check did not finish')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    walls=G.make_walls()
    signature=hashlib.sha256(b''.join(w.bez.tobytes() for w in walls)+
        str((DOMAIN,CELLS,RASTER,TOL)).encode()).hexdigest()
    f=C.Field3D('bspline',DOMAIN,CELLS,1)
    if (OUT/'field.npz').exists() and json.loads((OUT/'fit.json').read_text())['signature']==signature:
        f.set_W(np.load(OUT/'field.npz')['W'])
    else:
        print('Fitting bicubic workspace SDF (no robot footprint)...',flush=True)
        f=fit(walls)
        np.savez_compressed(OUT/'field.npz',W=f.W)
        (OUT/'fit.json').write_text(json.dumps(dict(signature=signature)))
    # Degenerate point body reduces the existing certificate to phi on wall curves.
    point=SimpleNamespace(bez=np.zeros((1,4,2)),rho=0.)
    union=SimpleNamespace(bez=np.concatenate([w.bez for w in walls]))
    cert=C.Certifier(f,point,union)
    level,ok,regions,seconds=cert.certify(n_theta0=1,tol=TOL,verbose=True)
    if not ok:
        raise RuntimeError('Boundary level could not be certified')
    interior=certify_interior(f,level,refined_walls(walls))
    band,bmin=cert.certify_band(level)
    regular,depth,remaining=cert.certify_regular(level,max_depth=12)
    if not (band and regular):
        raise RuntimeError('Field boundary/regularity checks failed')
    metadata=dict(signature=signature,level_m=level,tolerance_m=TOL,
        tightest_level_lower_m=level-TOL,tightest_level_upper_m=level,
        boundary_ok=bool(ok),band_ok=bool(band),regular_ok=bool(regular),
        regular_depth=depth,regular_open_boxes=remaining,band_min=float(bmin),
        cells=[CELLS,CELLS],raster_m=RASTER,regions=regions,certificate_s=seconds,
        **interior,interpretation='Workspace sublevel enclosure; no robot inflation.',
        arithmetic='Float64; no outward rounding.',
        source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                       for name in ('maze_levelset_figure.py','maze_reference_geometry.py')})
    (OUT/'results.json').write_text(json.dumps(metadata,indent=2))
    print(json.dumps(metadata,indent=2),flush=True)
    x=np.linspace(-4.8,4.8,1921); y=np.linspace(-5.55,5.55,2221)
    xx,yy=np.meshgrid(x,y)
    q=np.column_stack([xx.ravel(),yy.ravel(),np.zeros(xx.size)])
    values=np.concatenate([f.eval(q[i:i+10000]) for i in range(0,len(q),10000)]).reshape(xx.shape)
    np.savez_compressed(OUT/'displayed_field.npz',x=x,y=y,phi=values,level=level)
    plt.rcParams.update({'font.size':11,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10.2,6.4),sharex=True,sharey=True)
    for wall in walls:
        axes[0].add_patch(Polygon(wall.dense,facecolor='#b8bcc2',edgecolor='none'))
    axes[1].contourf(x,y,values,levels=[values.min()-1,level],colors=['#d4e1ef'])
    axes[1].contour(x,y,values,levels=[level],colors=['#143f70'],linewidths=1.4)
    axes[0].set_title('(a) Ground truth',pad=10)
    axes[1].set_title(r'(b) Fitted SDF: $\phi(x,y)=l$',pad=10)
    for ax in axes:
        ax.set(xlim=(-4.8,4.8),ylim=(-5.55,5.55),aspect='equal',xlabel='$x$ [m]')
        ax.spines[['top','right']].set_visible(False)
        ax.set_xticks(np.arange(-4,5,2));ax.set_yticks(np.arange(-4,5,2))
    axes[0].set_ylabel('$y$ [m]')
    fig.text(.75,.025,fr'$l={1000*level:.2f}$ mm; minimum enclosing level bracket $\leq {1000*TOL:.1f}$ mm',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1))
    fig.savefig(FIG/'maze_ground_truth_levelset.png',dpi=250,bbox_inches='tight')
    fig.savefig(FIG/'maze_ground_truth_levelset.pdf',bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    main()
