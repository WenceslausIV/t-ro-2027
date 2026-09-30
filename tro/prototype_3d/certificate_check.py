"""Independent cubic-field checks of the certificates, including negative multipliers.

These numerical checks detect implementation errors; they are not a safety proof.
Run: python prototype_3d/certificate_check.py
"""
import itertools
import importlib.util
import os
from pathlib import Path
for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import summed as S


class Cubic:
    def __init__(self, rng, dim):
        self.dim = dim
        self.v = rng.normal()
        self.g = rng.normal(size=dim)
        H = rng.normal(size=(dim, dim))
        self.H = (H + H.T) / 2
        T = rng.normal(size=(dim, dim, dim))
        self.T = sum(T.transpose(p) for p in itertools.permutations(range(3))) / 6
        self.M = np.linalg.norm(self.T)

    def eval(self, P, order=0):
        x = np.atleast_2d(P)[:, :self.dim]
        v = (self.v + x @ self.g + .5 * np.einsum('ni,ij,nj->n', x, self.H, x)
             + np.einsum('ijk,ni,nj,nk->n', self.T, x, x, x) / 6)
        if not order:
            return v
        g = self.g + x @ self.H + .5 * np.einsum('ijk,nj,nk->ni', self.T, x, x)
        H = self.H + np.einsum('ijk,nk->nij', self.T, x)
        gp = np.zeros((len(x), 3)); gp[:, :self.dim] = g
        Hp = np.zeros((len(x), 3, 3)); Hp[:, :self.dim, :self.dim] = H
        return (v, gp) if order == 1 else (v, gp, Hp)


class Bound:
    def __init__(self, M):
        self.M = M

    def eval(self, P):
        return np.full(len(P), self.M)


def check(dim):
    rng = np.random.default_rng(192 + dim)
    fa, fb = Cubic(rng, dim), Cubic(rng, dim)
    n, s = 24, .18
    r = s * np.sqrt(dim) / 2
    c = rng.normal(size=(n, dim)) * .3
    R, _ = np.linalg.qr(rng.normal(size=(dim, dim)))
    R[:, 0] *= np.linalg.det(R)
    yc = c @ R.T + rng.normal(size=dim) * .2
    va, ga, Ha = fa.eval(c, 2)
    D = dict(c=c, yc=yc, vA=va, gA=ga[:, :dim], HA=Ha[:, :dim, :dim], m3A=np.full(n, fa.M))
    val, g, H, m3 = S._layer(D, s, r, R, fb, 0., Bound(fb.M), dim, 5.)
    m = 3 if dim == 2 else 6
    W = np.zeros((m, dim, dim))
    if dim == 2:
        W[2] = [[0., -1.], [1., 0.]]
    else:
        for j in range(3):
            W[3 + j] = np.cross(np.eye(3)[j], np.eye(3)).T
    J = np.zeros((n, m, dim)); J[:, :dim] = np.eye(dim)
    J += np.einsum('mab,nb->nma', W, yc)
    Sb = np.c_[np.ones((n, dim)), np.repeat((np.linalg.norm(yc, axis=1) + r)[:, None], m - dim, axis=1)]
    baseline = Path(__file__).resolve().parent.parent / 'results/certificate_ablation/baseline/summed.py'
    if baseline.exists():
        spec = importlib.util.spec_from_file_location('summed_baseline', baseline)
        old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
        body = dict(PC=c, side=s, r=r, dim=dim, **{k: D[k] for k in ('vA', 'gA', 'HA', 'm3A')})
        args = (body, np.arange(n), R, yc, fb, 0., Bound(fb.M), J, W, Sb, 1e9, np.ones(m), np.ones(m))
        for row in ('vertex', 'bernstein'):
            old.ROW_MODE = row
            previous = old.rows(*args, depth=0, prune_rows=False)
            current = S.rows(*args, depth=0, prune_rows=False, row_mode=row, mult_mode='one')
            for a, b in zip(previous, current):
                np.testing.assert_array_equal(a, b)
        print(f'{dim}D unit multiplier exactly matches the saved original row builders.', flush=True)
    # Force the negative-projection case in addition to arbitrary orientations.
    assert np.any(S.multiplier_weights(D['gA'], g @ R, 'proj') < 0)
    count = 0
    for mult in S.MULT_MODES:
        w = S.multiplier_weights(D['gA'], g @ R, mult)
        cv = S.weighted_values(val, D, g, R, s, r, dim, mult)
        for mode, rowfn in (('vertex', S._vertex_rows), ('bernstein', S._velocity_rows)):
            A, T, C = rowfn(cv, g, H, m3, J, W, Sb, s, r, R, dim, 5.)
            K = 2 ** dim if mode == 'vertex' else 3 ** dim
            A, T, C = A.reshape(n, K, m), T.reshape(n, K, m), C.reshape(n, K)
            worst = np.inf
            for j in range(n):
                e = np.r_[rng.uniform(-s / 2, s / 2, (96, dim)),
                          np.array(list(itertools.product((-s / 2, s / 2), repeat=dim)))]
                x, y = c[j] + e, yc[j] + e @ R.T
                u = rng.uniform(-1, 1, (len(x), m))
                yd = u[:, :dim] + np.einsum('pm,mab,pb->pa', u, W, y)
                vb, gb = fb.eval(y, 1)
                actual = np.einsum('pi,pi->p', gb[:, :dim], yd) + 5 * (vb + w[j] * fa.eval(x))
                lower = (u @ A[j].T + np.abs(u) @ T[j].T + C[j]).min(axis=1)
                worst = min(worst, float((actual - lower).min()))
                count += len(x)
            assert worst >= -1e-10, (dim, mode, mult, worst)
            print(f'{dim}D {mode}/{mult}: min(actual - lower)={worst:.6g}', flush=True)
    return count


def negative_remainder_regression():
    # phi_A=x+x^3, phi_B=x, centered at zero. Projection gives w approximately -1.
    # Their linear models cancel, but the TRUE sum is negative for x>0. Using
    # w*M_A instead of abs(w)*M_A would falsely certify nonnegativity at rest.
    dim, side = 3, .2
    r = side * np.sqrt(dim) / 2
    D = dict(vA=np.zeros(1), gA=np.array([[1., 0., 0.]]),
             HA=np.zeros((1, dim, dim)), m3A=np.array([6.]))
    g = D['gA'].copy()
    base = 2 * side * S.TABLES[dim][0][None, :, 0] - r ** 3
    cv = S.weighted_values(base, D, g, np.eye(dim), side, r, dim, 'proj')
    w = S.multiplier_weights(D['gA'], g, 'proj')[0]
    x = side / 2
    actual = x + w * (x + x ** 3)
    assert cv.min() <= actual < 0
    incorrect = cv.min() + (abs(w) - w) * r ** 3
    assert incorrect > actual
    print('PASS: negative-multiplier cubic regression rejects the signed-remainder error.')


if __name__ == '__main__':
    negative_remainder_regression()
    print(f'PASS: {check(2) + check(3)} point/input comparisons (not a continuum proof).')
