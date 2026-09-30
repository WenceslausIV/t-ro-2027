"""Independent polygon intersection checks for recorded simulation poses."""
import numpy as np


def polygons_intersect(p, q):
    """Check edge intersections and containment, including crossings without inside vertices."""
    from matplotlib.path import Path
    if np.any(p.max(0) < q.min(0)) or np.any(q.max(0) < p.min(0)):
        return False
    a, b, c, d = p, np.roll(p, -1, axis=0), q, np.roll(q, -1, axis=0)
    candidates = ((np.minimum(a, b)[:, None] <= np.maximum(c, d)[None] + 1e-12) &
                  (np.maximum(a, b)[:, None] >= np.minimum(c, d)[None] - 1e-12)).all(axis=2)
    i, j = np.nonzero(candidates)

    def orient(a, b, c):
        v, w = b - a, c - a
        return v[:, 0] * w[:, 1] - v[:, 1] * w[:, 0]

    crossing = ((orient(a[i], b[i], c[j]) * orient(a[i], b[i], d[j]) <= 1e-24) &
                (orient(c[j], d[j], a[i]) * orient(c[j], d[j], b[i]) <= 1e-24))
    return bool(crossing.any() or Path(p).contains_point(q[0]) or Path(q).contains_point(p[0]))


def audit_trajectory(ground_truth, trajectory):
    """Check recorded poses only; this does not certify intersample collision avoidance."""
    from itertools import combinations
    contacts = []
    for k, poses in enumerate(trajectory):
        world = []
        for g, x in zip(ground_truth, poses):
            c, s = np.cos(x[2]), np.sin(x[2])
            world.append(g @ np.array([[c, s], [-s, c]]) + x[:2])
        for i, j in combinations(range(len(world)), 2):
            if polygons_intersect(world[i], world[j]):
                contacts.append([k, i, j])
    return dict(checked_poses=len(trajectory), intersection_count=len(contacts), contacts=contacts)
