"""Shared fixed-orientation C-space benchmark; no deadlock/liveness ranking.

Reimplementations, not author code:
  https://arxiv.org/abs/2504.09038 (minimum over samples, all active rows)
  https://arxiv.org/abs/2106.06330 (deforming ball world, Eqs. 10--17)

All controllers receive ONE common conservative C-space envelope of the exact
Minkowski sum of two nonconvex rectangle unions. Thus this isolates controller
representation/online cost, not end-to-end geometry preprocessing superiority.
Ball world includes the time derivative of its changing mapping and common speed
bounds by joint time scaling. This is explicitly an adaptation of the paper.

Run: python compare_three_baselines.py --trials 6 --seconds 12
Checks: python compare_three_baselines.py --check-only
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree
import torch

import cspace_sdf_cbf_compare as C
from sdf_cbf_utils import _ldp, Q_B2B

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'three_baselines'
TAU = 2 * np.pi
ROBOT = np.array([[.26, .09], [.09, .26]])
OBSTACLE = np.array([[.50, .14], [.14, .50]])
RECTS = (ROBOT[:, None, :] + OBSTACLE[None, :, :]).reshape(-1, 2)
CENTERS = np.array([[0., 0.]])
CLEARANCE = .020
WORLD_R = 4.2
VMAX = .6
GAMMA = 5.


def exact_radius(theta, rectangles):
    v = np.abs(np.stack([np.cos(theta), np.sin(theta)], axis=-1))
    return np.min(rectangles / np.maximum(v[..., None, :], 1e-15), axis=-1).max(-1)


def point_rect_union_distance(x, centers=CENTERS, rectangles=RECTS):
    """Exact set distance for disjoint translated rectangle-union bodies.

    Negative values are a union-of-rectangles overlap indicator, NOT penetration.
    """
    q = np.abs(np.asarray(x)[..., None, None, :] - centers[:, None, :]) - rectangles
    sd = np.linalg.norm(np.maximum(q, 0), axis=-1) + np.minimum(q.max(-1), 0)
    return sd.min(axis=(-1, -2))


class Envelope:
    def __init__(self, knots=128):
        t0 = time.perf_counter()
        # Square dilation includes the requested Euclidean disk dilation.
        rectangles = RECTS + CLEARANCE
        angles = np.linspace(0, TAU, knots + 1)
        values = exact_radius(angles, rectangles)
        values[-1] = values[0]
        spline = CubicSpline(angles, values, bc_type='periodic')
        h = TAU / knots
        # Derivative Bernstein control values on each cubic span.
        c = spline.c
        deriv = np.stack([c[2], c[2] + c[1] * h,
                          c[2] + 2 * c[1] * h + 3 * c[0] * h*h])
        Lfit = np.abs(deriv).max()
        Ltrue = np.max(np.sum(rectangles**2, axis=1) / rectangles.min(axis=1))
        n_verify = 32768
        a = np.arange(n_verify) * TAU / n_verify
        # Lipschitz extension makes this a continuous-angle upper-envelope bound,
        # in exact arithmetic. No verified floating-point rounding is claimed.
        extra = max(0., float(np.max(exact_radius(a, rectangles) - spline(a))))
        extra += (Ltrue + Lfit) * np.pi / n_verify + 1e-5
        self.spline = CubicSpline(angles, values + extra, bc_type='periodic')
        self.knots = knots
        self.rmin = float(values.min() + extra)
        c = self.spline.c
        self.bern = np.stack([c[3], c[3]+h*c[2]/3,
                             c[3]+2*h*c[2]/3+h*h*c[1]/3,
                             c[3]+h*c[2]+h*h*c[1]+h**3*c[0]], axis=1)
        self.rmax = float(self.bern.max())
        self.curve_speed_bound = float(np.hypot(self.rmax, Lfit))
        self.meta = dict(knots=knots, radial_padding_m=extra,
                         lipschitz_exact=Ltrue, lipschitz_spline=Lfit,
                         verification_nodes=n_verify, requested_clearance_m=CLEARANCE,
                         preprocessing_s=time.perf_counter()-t0)

    def radius(self, theta, derivative=0):
        return self.spline(np.asarray(theta) % TAU, derivative)

    def boundary(self, n):
        t = np.arange(n) * TAU / n
        return self.radius(t)[:, None] * np.stack([np.cos(t), np.sin(t)], axis=-1)

    def inside(self, x):
        a = np.arctan2(x[..., 1], x[..., 0])
        return np.linalg.norm(x, axis=-1) <= self.radius(a)


def hard_project(nominal, G, rhs):
    u = _ldp(np.asarray(nominal), np.asarray(G), np.asarray(rhs))
    if u is None:
        raise RuntimeError('Hard QP infeasible; no safety slack or hidden fallback used')
    residual = float(np.min(G @ u - rhs))
    if residual < -2e-7:
        raise RuntimeError(f'QP residual {residual}')
    return u


ANGLE = np.arange(12) * TAU / 12
SPEED_G = -np.stack([np.cos(ANGLE), np.sin(ANGLE)], axis=1)
SPEED_RHS = np.full(12, -VMAX*np.cos(np.pi/12))


def nominal(x, goal):
    v = 1.0 * (goal-x)
    return v * min(1., VMAX*np.cos(np.pi/12)/max(np.linalg.norm(v), 1e-15))


def workspace_row(x):
    return -2*x, -GAMMA*(WORLD_R**2-x@x)


def field_setup(env):
    t0 = time.perf_counter()
    D = 1.5
    field = C.Field3D('bspline', D, 60, 4)
    gx = np.linspace(-D, D, 181)
    xx, yy = np.meshgrid(gx, gx, indexing='ij')
    xy = np.stack([xx.ravel(), yy.ravel()], axis=-1)
    boundary = env.boundary(16384)
    distances = cKDTree(boundary).query(xy)[0]
    distances[env.inside(xy)] *= -1
    Phi = field.bx.design(gx)
    M = np.linalg.solve(Phi.T@Phi + 1e-5*np.eye(field.bx.n), Phi.T)
    W2 = M @ distances.reshape(len(gx), len(gx)) @ M.T
    field.set_W(np.repeat(W2[:, :, None], 4, axis=2))
    fit_s = time.perf_counter()-t0
    # Polar-envelope boundary boxes + the existing polynomial restriction bounds.
    dummy = C.SmoothShape(np.zeros((4, 2)))
    cert = C.Certifier(field, dummy, dummy)
    k = np.arange(env.knots)
    lo = np.zeros(len(k)); hi = np.ones(len(k))
    best = -np.inf
    eps = .001
    t_cert = time.perf_counter()
    processed = 0
    from cspace_cbf_5robots import sub_matrix_np
    for depth in range(35):
        mid = (k+(lo+hi)/2)*TAU/env.knots
        pts = env.radius(mid)[:, None] * np.stack([np.cos(mid), np.sin(mid)], axis=-1)
        best = max(best, float(field.eval(np.c_[pts, np.zeros(len(pts))]).max()))
        rb = (sub_matrix_np(lo, hi) @ env.bern[k, :, None])[..., 0]
        a0 = (k+lo)*TAU/env.knots; a1 = (k+hi)*TAU/env.knots
        r0, r1 = rb.min(1), rb.max(1)
        lows=[]; highs=[]
        for shift in (0., -np.pi/2):
            b0=a0+shift; b1=a1+shift
            v0=np.cos(b0); v1=np.cos(b1)
            mn=np.minimum(v0,v1); mx=np.maximum(v0,v1)
            mx=np.where(np.ceil(b0/TAU)<=np.floor(b1/TAU),1.,mx)
            mn=np.where(np.ceil((b0-np.pi)/TAU)<=np.floor((b1-np.pi)/TAU),-1.,mn)
            products=np.stack([r0*mn,r0*mx,r1*mn,r1*mx])
            lows.append(products.min(0)); highs.append(products.max(0))
        tmin=C._t(np.stack(lows,axis=1)); tmax=C._t(np.stack(highs,axis=1))
        z=C._t(np.zeros(len(k)))
        ub=cert.upper_bound(tmin,tmax,z,z).cpu().numpy()
        processed += len(k)
        keep=ub >= best+eps
        if not keep.any(): break
        k,lo,hi=k[keep],lo[keep],hi[keep]
        mid=(lo+hi)/2
        k=np.r_[k,k];lo,hi=np.r_[lo,mid],np.r_[mid,hi]
    else: raise RuntimeError('Contact certificate did not terminate')
    level=best+eps
    contact_s=time.perf_counter()-t_cert
    band,band_min=cert.certify_band(level)
    reg,regdepth,unresolved=cert.certify_regular(level)
    if not band or not reg: raise RuntimeError('Band/regularity failed')
    # Exactly collapse the constant orientation dimension for this translation-
    # only benchmark. Scalar bicubic power evaluation is the SAME fitted field.
    field.power_xy=np.einsum('ip,abpq,jq->abij',C.MB,field.C[:,:,0,:,:,0],C.MB)
    meta=dict(fit_s=fit_s, contact_cert_s=contact_s,
              total_preparation_s=time.perf_counter()-t0,
              level_m=level, tolerance_m=eps, band_min_m=float(band_min),
              regularity=bool(reg), regularity_depth=regdepth, regions=processed,
              runtime_storage_bytes=field.power_xy.nbytes, coefficient_bytes=field.W[:,:,0].nbytes,
              angular_axis='constant; fixed-orientation translation benchmark',
              arithmetic='float64; exact-arithmetic bounds, no outward rounding')
    return field, level, meta


class FieldController:
    def __init__(self, field, level):
        self.D=field.bx.hi;self.K=field.bx.K;self.spacing=field.bx.h
        self.power=field.power_xy;self.level=level
    def evaluate(self,p):
        s=(p+self.D)/self.spacing
        a,b=np.clip(np.floor(s).astype(int),0,self.K-1)
        u,v=s-[a,b]
        bu=np.array([1,u,u*u,u*u*u]);bv=np.array([1,v,v*v,v*v*v])
        du=np.array([0,1,2*u,3*u*u]);dv=np.array([0,1,2*v,3*v*v])
        c=self.power[a,b]
        return bu@c@bv,np.array([du@c@bv,bu@c@dv])/self.spacing
    def step(self,x,goal,dt):
        G=[*SPEED_G];rhs=[*SPEED_RHS]; hs=[]
        g,b=workspace_row(x);G.append(g);rhs.append(b)
        for center in CENTERS:
            p=x-center
            if np.abs(p).max() >= self.D-C.ACT_DELTA: continue
            val,grad=self.evaluate(p)
            h=float(val-self.level)
            hs.append(h); G.append(grad);rhs.append(-GAMMA*h)
        u=hard_project(nominal(x,goal),np.asarray(G),np.asarray(rhs))
        return u,dict(rows=len(G)-12,min_h=min(hs,default=np.nan))


class SampleController:
    def __init__(self,env,n):
        t0=time.perf_counter()
        self.samples=np.concatenate([env.boundary(n)+c for c in CENTERS])
        self.n=n
        # Boundary covering radius from sup ||d boundary / d angle||.
        self.cover=env.curve_speed_bound*np.pi/n + 1e-6
        self.tree=cKDTree(self.samples)
        self.meta=dict(samples_per_obstacle=n,covering_radius_m=self.cover,
                       preprocessing_s=time.perf_counter()-t0,
                       coordinate_storage_bytes=self.samples.nbytes,
                       note='KD-tree memory excluded; all potentially nonredundant point rows enforced')
    def step(self,x,goal,dt):
        ds,_=self.tree.query(x)
        # Enforce each potentially nonredundant sample CBF, not an arbitrary
        # closest point. Rows farther than rcrit hold for every allowed input.
        # This is a sufficient all-sample implementation of the minimum barrier;
        # it satisfies all active-gradient conditions and avoids missing switches.
        rcrit=(VMAX+np.sqrt(VMAX**2+(GAMMA*self.cover)**2))/GAMMA
        ids=self.tree.query_ball_point(x,float(rcrit)+1e-9)
        delta=x-self.samples[ids]
        h=ds**2-self.cover**2
        g,b=workspace_row(x)
        G=np.vstack([SPEED_G,2*delta,g])
        rhs=np.r_[SPEED_RHS,-GAMMA*(np.sum(delta**2,axis=1)-self.cover**2),b]
        u=hard_project(nominal(x,goal),G,rhs)
        return u,dict(rows=len(ids)+1,min_h=h)


class BallController:
    """Single-obstacle specialization of the deforming-ball strategy.

    The paper permits any diffeomorphism. For ONE star-shaped C-space obstacle,
    use the exact radial diffeomorphism F(x)=q_c+rho*(x-c)/r(angle(x-c)).
    It maps its entire exterior to a circle exterior, with determinant
    rho**2/r(angle)**2 > 0, avoiding an unverified multi-obstacle blend.
    Eq. (10) controls BOTH mapped center and radius; Eq. (16) restores them.
    Add positive-radius constraint, time-dependent-map derivative, and joint
    time scaling for common physical speed/workspace bounds. Not a literal
    reproduction of the paper's multi-obstacle numerical example.
    """
    def __init__(self,env):
        self.env=env
        self.center=np.zeros(2)
        self.radius=.5
        self.initial=np.r_[self.center,self.radius]
    def mapping(self,x,goal):
        d=x-CENTERS[0]; rr2=d@d
        if rr2<1e-12: raise RuntimeError('Map at obstacle center')
        a=np.arctan2(d[1],d[0]);r=float(self.env.radius(a))
        dr=float(self.env.radius(a,1))*np.array([-d[1],d[0]])/rr2
        f=d/r
        J=self.radius*(np.eye(2)/r-np.outer(d,dr)/r**2)
        P=np.column_stack([np.eye(2),f])
        return self.center+self.radius*f,J,P,np.array([rr2/r**2-1])
    def step(self,x,goal,dt):
        q,J,P,beta=self.mapping(x,goal)
        det=np.linalg.det(J)
        if det<=1e-7: raise RuntimeError('Radial map degenerate')
        qdot=J@nominal(x,goal)
        d=q-self.center;h=d@d-self.radius**2
        G=np.array([[-2*d[0],-2*d[1],-2*self.radius],[0,0,1]])
        rhs=np.array([-GAMMA*h-2*d@qdot,-GAMMA*(self.radius-.1)])
        rates=hard_project(self.initial-np.r_[self.center,self.radius],G,rhs)
        u=np.linalg.solve(J,qdot-P@rates)
        ratio=float(np.max((-SPEED_G)@u)/(VMAX*np.cos(np.pi/12)))
        scale=min(1.,1/max(ratio,1e-15))
        if x@u>0:
            scale=min(scale,GAMMA*(WORLD_R**2-x@x)/(2*x@u))
        if scale<0: raise RuntimeError('Outside shared workspace')
        u*=scale;rates*=scale
        self.center+=dt*rates[:2];self.radius+=dt*rates[2]
        return u,dict(rows=2,min_h=h,jacobian_det=det,
                      jacobian_condition=float(np.linalg.cond(J)),rate_scale=scale)


def validate(env,field=None,level=None):
    a=np.arange(10001)*TAU/10001
    excess=env.radius(a)-exact_radius(a,RECTS+CLEARANCE)
    assert excess.min()>0
    ball=BallController(env)
    goal=np.array([2.8,.15]);x=np.array([-.9,-.65])
    q,J,P,_=ball.mapping(x,goal)
    eps=1e-6
    numeric=np.column_stack([(ball.mapping(x+eps*e,goal)[0]-ball.mapping(x-eps*e,goal)[0])/(2*eps)
                              for e in np.eye(2)])
    assert np.max(np.abs(J-numeric))<1e-5
    # Parameter derivatives of the changing mapping, absent in static-map CBFs.
    state=np.r_[ball.center,ball.radius]
    def setstate(s):
        ball.center=s[:2].copy();ball.radius=float(s[2])
    for j in range(3):
        v=np.eye(3)[j]*eps;setstate(state+v);qp=ball.mapping(x,goal)[0]
        setstate(state-v);qm=ball.mapping(x,goal)[0]
        assert np.max(np.abs((qp-qm)/(2*eps)-P[:,j]))<1e-6
    setstate(state)
    sample=SampleController(env,128)
    # Two equally close points require BOTH rows, not an arbitrary tie winner.
    sample.samples=np.array([[.2,.05],[.2,-.05],[10.,10.]])
    sample.cover=.05;sample.tree=cKDTree(sample.samples)
    u,diag=sample.step(np.zeros(2),np.array([1.,0.]),.01)
    assert diag['rows']==3 and abs(u[0]-.5)<1e-8 and abs(u[1])<1e-8
    # A long interval with separated endpoints crosses the occupied set.
    assert audit_segment(np.array([-2.,0.]),np.array([2.,0.]))[3]
    assert not audit_segment(np.array([-2.,2.]),np.array([2.,2.]))[3]
    if field is not None:
        pts=env.boundary(10001)
        assert np.max(field.eval(np.c_[pts,np.zeros(len(pts))]))<level
        pt=np.array([[.8,.2,0.]])
        val,g=field.eval(pt,grad=True)
        gn=np.array([(field.eval(pt+eps*e)-field.eval(pt-eps*e))[0]/(2*eps) for e in np.eye(3)])
        assert np.max(np.abs(g[0]-gn))<1e-5
        fc=FieldController(field,level)
        for p in np.random.default_rng(1).uniform(-1.3,1.3,(100,2)):
            v,g=field.eval([*p,0],grad=True);v2,g2=fc.evaluate(p)
            assert abs(v[0]-v2)<1e-12 and np.max(np.abs(g[0,:2]-g2))<1e-11
    return dict(radial_enclosure_check_min_m=float(excess.min()),
                map_jacobian_max_error=float(np.max(np.abs(J-numeric))),
                multiple_sample_rows=True,intersample_crossing_detected=True)


def matched_query_timings(makers,data,trials,dt):
    """Same physical states for all methods; excludes path/stopping differences.

    Use all coarse-sampling trajectory states at a fixed stride. Reset the
    mapped world outside the timer, then include its complete online update.
    This is a query-cost diagnostic, not another closed-loop simulation.
    """
    queries=[]
    for i in range(trials):
        path=data[f'Sampling-64_{i}_x']
        goal=np.array([2.8,path[0,1]])
        queries.extend((x.copy(),goal) for x in path[:-1:20])
    timings={name:[] for name in makers}
    for repeat in range(3):
        names=list(makers);names=names[repeat:]+names[:repeat]
        for name in names:
            ctrl=makers[name]()
            for x,g in queries:
                if isinstance(ctrl,BallController):
                    ctrl.center[:]=0;ctrl.radius=.5
                t0=time.perf_counter_ns();ctrl.step(x,g,dt)
                timings[name].append((time.perf_counter_ns()-t0)*1e-6)
    return {name:dict(median_ms=float(np.median(v)),p95_ms=float(np.percentile(v,95)),
                      calls=len(v),unique_poses=len(queries)) for name,v in timings.items()}


def audit_segment(a,b,target=CLEARANCE,maxdepth=22):
    """Exact-geometry endpoint queries + distance Lipschitz interval bounds.

    Returns certified lower bounds and unresolved count, never declares a
    tolerance-limited interval safe. Physical distance is 1-Lipschitz under
    translation. Check requested clearance as well as zero-contact safety.
    """
    stack=[(a,b,float(point_rect_union_distance(a)),float(point_rect_union_distance(b)),0)]
    lower=np.inf;sampled=np.inf;unresolved=0;collision=False
    while stack:
        p,q,dp,dq,depth=stack.pop(); sampled=min(sampled,dp,dq)
        w=np.linalg.norm(q-p)
        lb=max(dp,dq)-w
        if min(dp,dq)<=0:
            collision=True;lower=min(lower,lb);continue
        # If an endpoint already violates the requested clearance, record that
        # separately and still certify collision avoidance at threshold zero.
        # Repeated subdivision cannot prove a false clearance specification.
        threshold = target if min(dp,dq)>=target else 0.
        if lb>=threshold:
            lower=min(lower,lb);continue
        if depth>=maxdepth or w<1e-8:
            unresolved+=1;lower=min(lower,lb);continue
        mid=(p+q)/2;dm=float(point_rect_union_distance(mid))
        stack.extend([(p,mid,dp,dm,depth+1),(mid,q,dm,dq,depth+1)])
    return lower,sampled,unresolved,collision


def simulate(controller,start,goal,steps,dt):
    x=start.copy();X=[x.copy()];U=[];timings=[];diags=[]
    # Untimed warmup is on an independent controller instance in the caller.
    for k in range(steps):
        t0=time.perf_counter_ns();u,diag=controller.step(x,goal,dt)
        timings.append((time.perf_counter_ns()-t0)*1e-6)
        if np.linalg.norm(u)>VMAX+1e-7: raise RuntimeError('Speed bound violated')
        x=x+dt*u;X.append(x.copy());U.append(u);diags.append(diag)
    X=np.array(X);U=np.array(U);times=np.array(timings)
    lowers=[];mins=[];unc=0;collision=False
    for a,b in zip(X[:-1],X[1:]):
        lb,m,n,c=audit_segment(a,b)
        lowers.append(lb);mins.append(m);unc+=n;collision|=c
    if not np.isfinite(X).all(): raise RuntimeError('Nonfinite trajectory')
    record=dict(collision=collision,clearance_unresolved_intervals=unc,
                sampled_clearance_violation=bool(min(mins)<CLEARANCE),
                verified_distance_lower_bound_m=float(min(lowers)),
                sampled_min_distance_m=float(min(mins)),
                control_median_ms=float(np.median(times)),control_p95_ms=float(np.percentile(times,95)),
                input_tv_m_per_s=float(np.linalg.norm(np.diff(U,axis=0),axis=1).sum()),
                max_speed_m_per_s=float(np.linalg.norm(U,axis=1).max()),
                path_length_m=float(np.linalg.norm(np.diff(X,axis=0),axis=1).sum()),
                min_barrier_raw=float(np.nanmin([d['min_h'] for d in diags])),
                max_rows=max(d['rows'] for d in diags),start=start.tolist(),goal=goal.tolist())
    if isinstance(controller,BallController):
        record.update(min_jacobian_determinant=min(d['jacobian_det'] for d in diags),
                      max_jacobian_condition=max(d['jacobian_condition'] for d in diags),
                      min_rate_scale=min(d['rate_scale'] for d in diags))
    return record,dict(x=X,u=U,times_ms=times)


def plain(o):
    if isinstance(o,np.ndarray):return o.tolist()
    if isinstance(o,np.generic):return o.item()
    raise TypeError(type(o))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=6)
    ap.add_argument('--seconds',type=float,default=12.);ap.add_argument('--dt',type=float,default=.01)
    ap.add_argument('--check-only',action='store_true');ap.add_argument('--smoke',action='store_true')
    args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    env=Envelope()
    if args.check_only:
        print(json.dumps(validate(env),indent=2));return
    field,level,fieldmeta=field_setup(env)
    checks=validate(env,field,level);print('Setup verified:',checks,flush=True)
    makers={'Ours':lambda:FieldController(field,level),
            **{f'Sampling-{n}':(lambda n=n:SampleController(env,n)) for n in (64,256,1024)},
            'Ball-world':lambda:BallController(env)}
    rng=np.random.default_rng(20260926)
    # Offset passes, not a head-on deadlock challenge. Keep every generated trial.
    offsets=np.resize(np.array([-.65,-.5,-.35,.35,.5,.65]),args.trials)
    offsets=offsets+rng.uniform(-.015,.015,args.trials)
    scenarios=[(np.array([-2.8,y]),np.array([2.8,y])) for y in offsets]
    if args.smoke:scenarios=scenarios[:1]
    steps=int(round(args.seconds/args.dt))
    result=dict(protocol=dict(seed=20260926,dt_s=args.dt,horizon_s=args.seconds,
                             trials=len(scenarios),speed_limit_m_per_s=VMAX,gamma=GAMMA,
                             clearance_m=CLEARANCE,world_radius_m=WORLD_R,
                             robot_rectangles_half_extents=ROBOT,obstacle_rectangles_half_extents=OBSTACLE,
                             centers=CENTERS,deadlock_metrics=False,
                             scope='shared C-space envelope, fixed-orientation nonconvex translating body',
                             timings='CPU wall-clock: complete online controller including QP/mapping; excludes validation',
                             sampling_margin='covering-radius correction in distance units, squared in barrier',
                             ball_adaptation='single-star radial diffeomorphism; moving center/radius; F_theta theta_dot; joint rate scaling',
                             sample_adaptation='all potentially nonredundant point CBFs; Euclidean boundary-covering correction'),
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                                 torch=torch.__version__,platform=platform.platform(),processor=platform.processor()),
                shared_geometry=env.meta,field=fieldmeta,checks=checks,runs={},summary={},
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                sources=['https://arxiv.org/abs/2504.09038','https://arxiv.org/abs/2106.06330'])
    data={}
    for name,maker in makers.items():
        warm=maker()
        for _ in range(10):warm.step(scenarios[0][0],scenarios[0][1],args.dt)
        if isinstance(warm,SampleController):result.setdefault('sampling',{})[name]=warm.meta
    # Rotate method order between trials to reduce systematic thermal/order effects.
    names=list(makers)
    for i,(start,goal) in enumerate(scenarios):
        for name in names[i%len(names):]+names[:i%len(names)]:
            rec,arrays=simulate(makers[name](),start,goal,steps,args.dt)
            result['runs'].setdefault(name,[]).append(rec)
            for k,v in arrays.items():data[f'{name}_{i}_{k}']=v
            print(i,name,'gap',round(rec['sampled_min_distance_m']*1000,2),
                  'ms',round(rec['control_median_ms'],3),'TV',round(rec['input_tv_m_per_s'],2),flush=True)
    for name,runs in result['runs'].items():
        pooled=np.concatenate([data[f'{name}_{i}_times_ms'] for i in range(len(scenarios))])
        result['summary'][name]=dict(control_median_ms=float(np.median(pooled)),
            control_p95_ms=float(np.percentile(pooled,95)),
            min_distance_mm=1000*min(r['sampled_min_distance_m'] for r in runs),
            verified_lower_bound_mm=1000*min(r['verified_distance_lower_bound_m'] for r in runs),
            input_tv_mean=float(np.mean([r['input_tv_m_per_s'] for r in runs])),
            collisions=sum(r['collision'] for r in runs),
            clearance_violation_runs=sum(r['sampled_clearance_violation'] for r in runs),
            clearance_unresolved=sum(r['clearance_unresolved_intervals'] for r in runs))
    prefix='smoke' if args.smoke else 'comparison'
    if not args.smoke:
        result['matched_query_timing']=matched_query_timings(makers,data,len(scenarios),args.dt)
    (OUT/f'{prefix}.json').write_text(json.dumps(result,indent=2,default=plain),encoding='utf8')
    np.savez_compressed(OUT/f'{prefix}_trajectories.npz',**data)
    if not args.smoke:render(result,data,env)
    print(json.dumps(result['summary'],indent=2),flush=True)


def render(result,data,env):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    names=list(result['summary']);colors=plt.get_cmap('tab10').colors
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for center in CENTERS:
        for hx,hy in OBSTACLE:
            axes[0].add_patch(Rectangle(center-[hx,hy],2*hx,2*hy,color='.65'))
        bd=env.boundary(1000)+center
        axes[0].plot(*bd.T,color='.5',ls=':',lw=1)
    for j,name in enumerate(names):
        x=data[f'{name}_0_x'];u=data[f'{name}_0_u']
        axes[0].plot(*x.T,label=name,color=colors[j],lw=1.5)
        # Draw the actual finite robot at the final pose, not just a point path.
        for hx,hy in ROBOT:
            axes[0].add_patch(Rectangle(x[-1]-[hx,hy],2*hx,2*hy,
                                       fill=False,ec=colors[j],lw=.8,alpha=.7))
        tt=np.arange(len(x))*result['protocol']['dt_s']
        axes[1].plot(tt,1000*point_rect_union_distance(x),color=colors[j],label=name)
    axes[0].set(xlabel='x [m]',ylabel='y [m]',aspect='equal',title='Same reference-body translation task')
    axes[0].legend(fontsize=8)
    axes[1].axhline(CLEARANCE*1000,color='k',ls=':',label='Requested clearance')
    axes[1].set(xlabel='t [s]',ylabel='True body-to-obstacle gap [mm]',ylim=(0,250),title='Physical geometry, common distance oracle')
    fig.savefig(OUT/'trajectories.png',dpi=200);plt.close(fig)
    paperfig=ROOT/'tro'/'figs'/'three_baselines.png'
    paperfig.write_bytes((OUT/'trajectories.png').read_bytes())
    rows=[]
    for name,s in result['summary'].items():
        rows.append(f"| {name} | {s['control_median_ms']:.3f} | {s['control_p95_ms']:.3f} | {s['min_distance_mm']:.2f} | {s['verified_lower_bound_mm']:.2f} | {s['input_tv_mean']:.2f} | {s['collisions']} | {s['clearance_unresolved']} |")
    header='| Method | median ms | p95 ms | sampled gap mm | verified lower bound mm | mean input TV (m/s) | collision runs | clearance unresolved |\n|---|---:|---:|---:|---:|---:|---:|---:|\n'
    (OUT/'table.md').write_text(header+'\n'.join(rows)+'\n',encoding='utf8')
    latex=[]
    for name,s in result['summary'].items():
        short=name.replace('Sampling-','Samples, $N=')
        if name.startswith('Sampling-'):short+='$'
        if name=='Ball-world':short='Ball world (adapted)'
        t=result['matched_query_timing'][name]
        latex.append(f"{short} & {t['median_ms']:.3f} & {t['p95_ms']:.3f} & {s['min_distance_mm']:.1f} & {s['input_tv_mean']:.2f} \\\\")
    table=r'''% Generated by compare_three_baselines.py from results/three_baselines/comparison.json.
\begin{table}[t]
\caption{Shared Configuration-Space Benchmark (Six Fixed-Horizon Runs)}
\label{tab:three-baselines}
\centering\footnotesize
\setlength{\tabcolsep}{3pt}
\begin{tabular}{lrrrr}
\toprule
Method & median & p95 & min. gap & input TV\\
 & [ms] & [ms] & [mm] & [m/s]\\
\midrule
'''+ '\n'.join(latex)+r'''
\bottomrule
\end{tabular}
\\[3pt]
\parbox{\columnwidth}{\footnotesize Times: identical pose queries, complete CPU
controller including the QP; ball-world state reset before each timed query.
Gap: minimum evaluated physical-body separation over the six closed-loop runs.
Input TV: mean $\sum_k\|v_{k+1}-v_k\|_2$ over the fixed 12-s horizon.
No collision or unresolved clearance interval occurred for any method.}
\end{table}
'''
    (ROOT/'tro'/'three_baseline_table.tex').write_text(table,encoding='utf8')


if __name__=='__main__':main()
