"""
Ground-truth outlines for the docking experiment (paper/main.tex, "Docking into a tapered receptacle").

    station_body()   : module with radiator panels and a tapered receptacle opening towards +x
    craft_body()     : vehicle with a long tapered nose along +x, panel arrays on booms, two nozzles

The receptacle is the seated nose offset outward by CLEARANCE everywhere: along the tapered walls,
in front of the nose tip, and between the two front faces; its corner radii follow those of the
nose. With CLEARANCE equal to the stopping distance of the certified barrier, the vehicle stops
fully seated with a thin, uniform gap. Both outlines are non-convex CCW polygons; certified
B-spline boundaries are fitted with cspace_cbf_5robots.certify_curve_encloses.
"""
import numpy as np

CLEARANCE = 0.032      # gap between the seated nose and the receptacle, all around [m]
MOUTH_X = 0.35         # station: x of the receptacle mouth (station front face)
NOSE_BASE_X = 0.30     # craft: x of the nose base (= front face of the hull)
NOSE_LEN = 0.70        # craft: nose length
NOSE_TIP_X = NOSE_BASE_X + NOSE_LEN
W_BASE, W_TIP = 0.16, 0.05                 # nose half-width at the base and at the tip
R_TIP = 0.04           # fillet radius of the nose-tip corners
R_BASE = 0.035         # fillet radius of the (concave) corners where the nose meets the hull
SEAT_DEPTH = NOSE_LEN - CLEARANCE          # tip depth past the mouth at full seating
SEAT_X = MOUTH_X + CLEARANCE + NOSE_BASE_X  # station-frame x of the craft origin at full seating (heading pi)


def fillet(P, r_convex=0.05, r_concave=0.015, n_arc=10, radius=None):
    """Round every corner of a CCW polygon: circular arcs of radius r_convex at convex corners and
    r_concave at concave ones (each limited to 45% of the adjacent edges); radius[k], if given and
    not None, overrides the radius at vertex k."""
    P = np.asarray(P, dtype=np.float64)
    out = []
    n = len(P)
    for k in range(n):
        u, v, w = P[k - 1], P[k], P[(k + 1) % n]
        d1, d2 = u - v, w - v
        l1, l2 = np.linalg.norm(d1), np.linalg.norm(d2)
        d1, d2 = d1 / l1, d2 / l2
        cross = (v - u)[0] * (w - v)[1] - (v - u)[1] * (w - v)[0]
        alpha = np.arccos(np.clip(d1 @ d2, -1.0, 1.0))
        if alpha > np.pi - 1e-6:                         # straight vertex
            out.append(v)
            continue
        r = r_convex if cross > 0 else r_concave
        if radius is not None and radius[k] is not None:
            r = radius[k]
        if r <= 0:
            out.append(v)
            continue
        t = min(r / np.tan(alpha / 2), 0.45 * min(l1, l2))
        r = t * np.tan(alpha / 2)
        bis = (d1 + d2) / np.linalg.norm(d1 + d2)
        c = v + bis * r / np.sin(alpha / 2)
        p1, p2 = v + d1 * t, v + d2 * t
        a1 = np.arctan2(*(p1 - c)[::-1])
        a2 = np.arctan2(*(p2 - c)[::-1])
        da = (a2 - a1 + np.pi) % (2 * np.pi) - np.pi     # short way round
        for s in np.linspace(0, 1, n_arc):
            out.append(c + r * np.array([np.cos(a1 + s * da), np.sin(a1 + s * da)]))
    return np.array(out)


def densify(P, ds=0.004):
    P = np.asarray(P, dtype=np.float64)
    out = []
    for a, b in zip(P, np.roll(P, -1, axis=0)):
        n = max(int(np.ceil(np.linalg.norm(b - a) / ds)), 1)
        out.append(a + (np.arange(n) / n)[:, None] * (b - a))
    return np.concatenate(out)


def _mirror(upper):
    return [(x, -y) for (x, y) in reversed(upper)]


def station_body():
    g = CLEARANCE
    slope = (W_BASE - W_TIP) / NOSE_LEN
    face = MOUTH_X + g                          # seated craft front face (station frame)
    x_tip = face - NOSE_LEN                     # seated nose tip
    wall = g / np.cos(np.arctan(slope))         # vertical offset giving a normal gap g
    w_mouth = W_BASE - slope * g + wall         # receptacle half-width at the mouth
    w_bot = W_TIP + wall
    upper = [(x_tip - g, w_bot), (MOUTH_X, w_mouth), (MOUTH_X, 0.30), (0.28, 0.42),
             (-0.54, 0.42), (-0.54, 0.52), (-0.38, 0.52), (-0.38, 0.98), (-0.80, 0.98), (-0.80, 0.52),
             (-0.64, 0.52), (-0.64, 0.42), (-0.90, 0.42)]
    # radii: pocket corners follow the nose tip (R_TIP + g), mouth corners the nose base (R_BASE - g)
    r_up = [R_TIP + g, max(R_BASE - g, 0.0)] + [None] * (len(upper) - 2)
    P = upper + _mirror(upper)
    return densify(fillet(P, radius=r_up + r_up[::-1]))


def craft_body():
    upper = [(NOSE_TIP_X, W_TIP), (NOSE_BASE_X, W_BASE), (NOSE_BASE_X, 0.20), (-0.04, 0.20),
             (-0.04, 0.30), (0.14, 0.30), (0.14, 0.72), (-0.30, 0.72), (-0.30, 0.30), (-0.12, 0.30),
             (-0.12, 0.20), (-0.40, 0.20), (-0.40, 0.16), (-0.52, 0.19), (-0.52, 0.03), (-0.40, 0.06)]
    r_up = [R_TIP, R_BASE] + [None] * (len(upper) - 2)
    P = upper + _mirror(upper)
    return densify(fillet(P, radius=r_up + r_up[::-1]))
