"""Vector reconstruction of the user-supplied winding maze (decorations omitted).

Image coordinates below are manually traced landmarks, not pixel-exact image extraction.
The interpolated piecewise cubic wall regions define the simulation geometry.
"""
import numpy as np
from matplotlib.path import Path
from scipy.spatial import cKDTree
import cspace_sdf_cbf_compare as C

SCALE = .03
ORIGIN = np.array([181.,238.])


def world(points):
    return (np.asarray(points,dtype=float)-ORIGIN)*[SCALE,-SCALE]


# The two long opposing boundaries follow the entrance, winding central passage,
# curled blind branches, and lower-right exit of the reference.
RIGHT_TRACE = [
    (97,65),(85,93),(76,124),(74,145),(79,161),(91,168),(105,166),
    (121,154),(138,131),(153,115),(173,106),(198,102),(219,104),(231,107),
    (246,88),(264,79),(279,79),(295,87),(307,101),(313,117),(309,131),
    (299,138),(289,136),(283,128),(284,120),(291,116),(282,105),(272,102),
    (261,106),(252,116),(266,130),(277,147),(283,169),(282,188),(275,203),
    (261,211),(242,216),(225,222),(214,230),(232,237),(251,244),(264,252),
    (274,249),(284,247),(291,241),(294,229),(285,226),(280,219),(283,210),
    (292,207),(303,212),(311,222),(315,236),(312,252),(305,265),(295,271),
    (284,274),(288,290),(285,305),(279,317),(289,328),(297,343),(301,357),
    (296,364),(286,363),(279,356),(275,342),(268,334),(257,332),(241,336),
    (222,336),(203,333),(185,325),(174,338),(156,348),(136,352),(120,352),
    (113,346),(112,337),(117,330),(132,329),(146,325),(162,313),(148,307),
    (128,303),(112,307),(101,317),(96,331),(96,347),(103,361),(115,370),
    (133,375),(158,377),(185,376),(212,374),(239,372),(260,376),(282,384),
    (303,388),(313,392),(314,411),
]
LEFT_TRACE = [
    (33,65),(43,86),(51,106),(52,125),(51,144),(54,164),(61,180),
    (73,188),(90,191),(109,190),(123,200),(132,213),(136,227),(132,238),
    (124,244),(108,241),(88,236),(73,236),(58,241),(47,253),(41,268),
    (41,282),(46,297),(53,303),(60,301),(64,295),(62,284),(62,275),
    (69,264),(79,259),(89,260),(106,265),(124,267),(139,262),(151,252),
    (157,238),(155,221),(148,203),(139,189),(126,179),(139,166),(152,146),
    (168,132),(189,125),(209,125),(229,133),(247,146),(261,163),
    (244,151),(229,143),(214,140),(201,143),(190,151),(185,162),(187,174),
    (196,183),(209,187),(218,183),(221,176),(217,168),(213,164),(222,164),
    (234,169),(245,181),(253,191),(235,195),(217,200),(203,207),(194,218),
    (190,231),(194,243),(207,253),(224,259),(242,266),(256,275),(264,287),
    (265,298),(257,307),(244,312),(230,313),(213,310),(195,302),(176,293),
    (155,285),(134,281),(113,281),(96,285),(84,293),(77,306),(73,322),
    (73,339),(77,357),(85,373),(98,385),(115,393),(137,398),(163,400),
    (189,400),(212,397),(231,395),(240,398),(244,403),(244,411),
]


class TracedWall:
    def __init__(self, trace, closure):
        pts = world(trace)
        delta = np.diff(pts,axis=0)
        lengths = np.linalg.norm(delta,axis=1)
        tangent = np.empty_like(pts)
        tangent[0],tangent[-1] = delta[0]/lengths[0],delta[-1]/lengths[-1]
        tangent[1:-1] = (pts[2:]-pts[:-2])/(lengths[:-1]+lengths[1:])[:,None]
        # Chord-length cubic Hermite interpolation preserves the hand-traced curls.
        bez = [np.array([a,a+h*ta/3,b-h*tb/3,b])
               for a,b,ta,tb,h in zip(pts[:-1],pts[1:],tangent[:-1],tangent[1:],lengths)]
        closing = np.vstack([pts[-1],world(closure),pts[0]])
        for a,b in zip(closing[:-1],closing[1:]):
            bez.append(np.array([a,(2*a+b)/3,(a+2*b)/3,b]))
        self.bez = np.array(bez)
        # Certified (exact-arithmetic) chord deviation bound for a cubic:
        # sup ||B''||/(8*n^2). Use <=0.1 mm and <=4 cm straight segments.
        pieces, errors = [],[]
        for span in self.bez:
            second = 6*np.diff(span,n=2,axis=0)
            m2 = np.linalg.norm(second,axis=1).max()
            n = max(4,int(np.ceil(np.sqrt(m2/(8e-4)))),
                    int(np.ceil(np.linalg.norm(np.diff(span,axis=0),axis=1).sum()/.04)))
            pieces.append(C.bern3(np.arange(n)/n)@span)
            errors.append(m2/(8*n*n))
        self.dense = self.outline = np.vstack(pieces)
        self.chord_error = float(max(errors))
        self.path = Path(self.dense)
        self.edge_a = self.dense
        self.edge_b = np.roll(self.dense,-1,axis=0)
        self.mid_tree = cKDTree((self.edge_a+self.edge_b)/2)
        self.vertex_tree = cKDTree(self.dense)
        self.max_half_edge = np.linalg.norm(self.edge_b-self.edge_a,axis=1).max()/2


class ContourWall(TracedWall):
    """Closed smooth wall from a closed world-frame contour (periodic chord-length Hermite cubics).
    Used for the fillers that close the maze's dead ends."""
    def __init__(self, contour, spacing=.04):
        c = np.asarray(contour,dtype=float)
        if np.linalg.norm(c[0]-c[-1]) < 1e-9:
            c = c[:-1]
        closed = np.vstack([c,c[:1]])
        arc = np.r_[0,np.linalg.norm(np.diff(closed,axis=0),axis=1).cumsum()]
        m = max(8,int(arc[-1]/spacing))
        q = np.linspace(0,arc[-1],m,endpoint=False)
        pts = np.column_stack([np.interp(q,arc,closed[:,j]) for j in range(2)])
        nxt, prv = np.roll(pts,-1,axis=0), np.roll(pts,1,axis=0)
        lengths = np.linalg.norm(nxt-pts,axis=1)
        tangent = (nxt-prv)/(lengths+np.roll(lengths,1))[:,None]
        self.bez = np.array([np.array([a,a+h*ta/3,b-h*tb/3,b])
                             for a,b,ta,tb,h in zip(pts,nxt,tangent,np.roll(tangent,-1,axis=0),lengths)])
        pieces, errors = [],[]
        for span in self.bez:
            second = 6*np.diff(span,n=2,axis=0)
            m2 = np.linalg.norm(second,axis=1).max()
            n = max(4,int(np.ceil(np.sqrt(m2/(8e-4)))),
                    int(np.ceil(np.linalg.norm(np.diff(span,axis=0),axis=1).sum()/.04)))
            pieces.append(C.bern3(np.arange(n)/n)@span)
            errors.append(m2/(8*n*n))
        self.dense = self.outline = np.vstack(pieces)
        self.chord_error = float(max(errors))
        self.path = Path(self.dense)
        self.edge_a = self.dense
        self.edge_b = np.roll(self.dense,-1,axis=0)
        self.mid_tree = cKDTree((self.edge_a+self.edge_b)/2)
        self.vertex_tree = cKDTree(self.dense)
        self.max_half_edge = np.linalg.norm(self.edge_b-self.edge_a,axis=1).max()/2


def make_walls():
    return [TracedWall(RIGHT_TRACE,[(329,411),(329,65)]),
            TracedWall(LEFT_TRACE,[(33,411)])]


def polygon_gap(pose, robot_gt, wall):
    """Exact polygon distance with conservative spatial pruning, less curve error."""
    p = robot_gt@C.rot(pose[2]).T+pose[:2]
    q = np.roll(p,-1,axis=0)
    rho = np.linalg.norm(robot_gt,axis=1).max()
    upper = wall.vertex_tree.query(p)[0].min()
    ids = wall.mid_tree.query_ball_point(pose[:2],rho+upper+wall.max_half_edge+1e-10)
    a,b = wall.edge_a[ids],wall.edge_b[ids]
    def cross(x,y):
        return x[...,0]*y[...,1]-x[...,1]*y[...,0]
    pq,ab = q-p,b-a
    # Distances of each endpoint to the opposite edge family.
    v = p[:,None]-a[None]
    t = np.clip(np.sum(v*ab[None],axis=2)/np.sum(ab*ab,axis=1)[None],0,1)
    da = np.linalg.norm(v-t[:,:,None]*ab[None],axis=2).min()
    v2 = a[None]-p[:,None]
    s = np.clip(np.sum(v2*pq[:,None],axis=2)/np.sum(pq*pq,axis=1)[:,None],0,1)
    db = np.linalg.norm(v2-s[:,:,None]*pq[:,None],axis=2).min()
    den = cross(pq[:,None],ab[None])
    good = abs(den)>1e-14
    u = np.divide(cross(v2,ab[None]),den,out=np.full_like(den,np.inf),where=good)
    v = np.divide(cross(v2,pq[:,None]),den,out=np.full_like(den,np.inf),where=good)
    hit = np.any((u>=0)&(u<=1)&(v>=0)&(v<=1))
    if hit or wall.path.contains_point(p[0]) or Path(p).contains_points(a).any():
        return -wall.chord_error
    return float(min(da,db)-wall.chord_error)


def physical_gaps(pose, robot_gt, walls):
    return np.array([polygon_gap(pose,robot_gt,w) for w in walls])


class FittedWall:
    """Nonuniform periodic cubic B-spline, exposed in cellwise Bernstein form."""
    def __init__(self, knots, points):
        from scipy.interpolate import CubicSpline, BSpline
        self.spline = CubicSpline(knots, points, bc_type='periodic')
        basis = [BSpline.from_power_basis(CubicSpline(knots,points[:,j],bc_type='periodic'),
                                         bc_type='periodic') for j in range(2)]
        self.ctrl, self.knots = np.column_stack([b.c for b in basis]), basis[0].t
        h = np.diff(knots)[:,None]
        ends = self.spline(knots)
        deriv = self.spline(knots, nu=1)
        self.bez = np.stack([ends[:-1],ends[:-1]+h*deriv[:-1]/3,
                             ends[1:]-h*deriv[1:]/3,ends[1:]],axis=1)
        self.span_lengths = np.diff(knots)
        self.dense = np.vstack([C.bern3(np.arange(n)/n)@b for b in self.bez
            for n in [max(3,int(np.ceil(np.linalg.norm(np.diff(b,axis=0),axis=1).sum()/.003)))]] )


def _offset_samples(wall, offset):
    """Round convex joins; intersect offset lines at concave polygon vertices."""
    p = wall.dense
    area = np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1))
    if area < 0:
        p = p[::-1]
    edges = np.roll(p,-1,axis=0)-p
    tang = edges/np.linalg.norm(edges,axis=1)[:,None]
    normals = np.column_stack([tang[:,1],-tang[:,0]])
    pieces = []
    for i,point in enumerate(p):
        u,v = tang[i-1],tang[i]
        n0,n1 = normals[i-1],normals[i]
        turn = np.arctan2(u[0]*v[1]-u[1]*v[0],u@v)
        if turn > .03:
            angle = np.arctan2(n0[1],n0[0])
            angles = angle+np.linspace(0,turn,max(2,int(np.ceil(turn/.12))+1))
            pieces.extend(point+offset*np.column_stack([np.cos(angles),np.sin(angles)]))
        else:
            pieces.append(point+offset*(n0+n1)/(1+n0@n1))
    points = np.array(pieces)
    points = points[np.linalg.norm(points-np.roll(points,1,axis=0),axis=1)>1e-8]
    points = np.vstack([points,points[:1]])
    # Preserve short corner arcs while adding samples along long straight spans.
    dense = []
    for a,b in zip(points[:-1],points[1:]):
        n=max(1,int(np.ceil(np.linalg.norm(b-a)/.025)))
        dense.extend(a+(b-a)*np.arange(n)[:,None]/n)
    closed = np.vstack([dense,dense[:1]])
    arc = np.r_[0,np.linalg.norm(np.diff(closed,axis=0),axis=1).cumsum()]
    return arc,closed


def boundary_distance_upper(points,wall):
    # Distance to a subset of polygon edges upper-bounds distance to its boundary.
    _,ids = wall.mid_tree.query(points,k=8)
    a=wall.edge_a[ids]; v=wall.edge_b[ids]-a; d=points[:,None]-a
    t=np.clip(np.sum(d*v,axis=2)/np.sum(v*v,axis=2),0,1)
    return np.linalg.norm(d-t[:,:,None]*v,axis=2).min(axis=1)


def fit_wall(wall):
    """Fit tightly, then check continuous enclosure and a 3-mm distance bound.

    Corner arcs get short nonuniform knot spans so a sharp corner does not force
    a large offset along the rest of the wall. Float64, without outward rounding.
    """
    from cspace_cbf_5robots import sub_matrix_np, winding_number
    left,right = sub_matrix_np(np.array([0.,.5]),np.array([.5,1.]))
    def box(b):
        lo,hi=b.min(axis=0)-1e-12,b.max(axis=0)+1e-12
        return np.array([lo,[hi[0],lo[1]],hi,[lo[0],hi[1]]])
    for offset in (.001,.0015,.002):
        arc,points=_offset_samples(wall,offset)
        shape=FittedWall(arc,points)
        ordered=[]; min_separation=np.inf; ok=True
        for b in shape.bez:
            stack=[(b,0)]
            while stack:
                part,depth=stack.pop()
                hull=box(part); center=hull.mean(axis=0)
                separation=polygon_gap(np.r_[center,0.],hull-center,wall)
                if separation>1e-8:
                    ordered.append(part[0]); min_separation=min(min_separation,separation)
                elif depth<12:
                    stack.extend([(right@part,depth+1),(left@part,depth+1)])
                else:
                    ok=False; break
            if not ok:
                break
        if ok:
            ok=winding_number(np.array(ordered),wall.dense[0])!=0
        print(f'Wall fit: {len(shape.ctrl)} controls, offset {offset:.4f} m, enclosure {ok}',flush=True)
        if not ok:
            continue
        # Every curve point is within speed/(2*n) <= 0.25 mm of a sample.
        # Include reference chord error to bound distance to the original curve.
        samples=[]
        for b in shape.bez:
            speed=3*np.linalg.norm(np.diff(b,axis=0),axis=1).max()
            n=max(1,int(np.ceil(speed/.0005)))
            samples.extend(C.bern3(np.linspace(0,1,n+1))@b)
        distances=boundary_distance_upper(np.array(samples),wall)
        upper=float(distances.max()+.00025+wall.chord_error)
        print(f'Continuous wall-fit distance upper bound: {upper*1000:.3f} mm',flush=True)
        if upper>.003:
            raise RuntimeError('Fitted boundary exceeds the 3-mm tightness requirement')
        shape.chord_error=wall.chord_error
        return shape,dict(controls=len(shape.ctrl),spans=len(shape.bez),offset_m=offset,
            enclosure_ok=True,hull_separation_lower_m=float(min_separation),
            boundary_distance_upper_m=upper,sampled_boundary_distance_max_m=float(distances.max()),
            sampled_boundary_distance_median_m=float(np.median(distances)),basis='nonuniform periodic cubic B-spline')
    raise RuntimeError('Could not certify a tight fitted wall enclosure')
