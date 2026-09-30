"""Analytic counterexample to inferring solid enclosure from boundary values.

phi is a C2 tensor-product cubic spline: a quadratic plus a tensor product of
centered cardinal cubic B-splines. The ground truth is the closed unit disk.
Sampling is a diagnostic only; the printed boundary bound is analytic.
"""
import json
from pathlib import Path
import numpy as np


def beta(t):
    a = np.abs(t)
    return np.where(a < 1, 2 / 3 - a * a + .5 * a ** 3,
                    np.where(a < 2, (2 - a) ** 3 / 6, 0.))


def phi(x, y):
    return x * x + y * y - 1 + 4.5 * beta(2 * x) * beta(2 * y)


def main():
    level = .2
    # On the unit circle, at least one |coordinate| >= 1/sqrt(2).
    # beta decreases with |t|, and beta <= 2/3 everywhere.
    upper = 4.5 * ((2 - np.sqrt(2)) ** 3 / 6) * (2 / 3)
    assert upper < level < phi(0., 0.)
    # phi >= x^2+y^2-1, so its level sublevel is compact inside [-2,2]^2.
    t = np.linspace(0, 2 * np.pi, 10001)
    result = dict(level=level, analytic_boundary_upper=float(upper),
                  diagnostic_sample_max=float(phi(np.cos(t), np.sin(t)).max()),
                  center_value=float(phi(0., 0.)),
                  sublevel_radius_upper=float(np.sqrt(1 + level)),
                  conclusion='Strict boundary inclusion does not imply solid enclosure, even for a C2 cubic spline with a compact sublevel set.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
