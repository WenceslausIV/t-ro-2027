"""
Goal marker shared by the Franka and maze figures: an exit-sign green ball with a white G wrapped onto its
surface (the letter lies on the sphere, so it bends with the curvature). Shading as in franka_gif.shade:
k = .45 + .55 |n . LIGHT| under the Franka figures' light and camera.

    ball_triangles(center, radius, elev, azim)  -> triangles and RGBA colors for a Poly3DCollection
    ball_image(elev, azim)                      -> RGBA image of the same ball as seen by that camera
"""
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path
from matplotlib.textpath import TextPath

GREEN, WHITE = np.array([0., .65, .32]), np.array([1., 1., 1.])
LIGHT = np.array([.4, -.5, .77])                      # franka_gif.LIGHT
G_HEIGHT = 1.0                                        # angular height of the letter on the sphere [rad]


def camera(elev, azim):
    """Unit vectors (toward viewer, screen right, screen up) of a matplotlib 3D view."""
    e, a = np.radians(elev), np.radians(azim)
    d = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    right = np.cross(-d, [0., 0., 1.])
    right /= np.linalg.norm(right)
    return d, right, np.cross(right, -d)


def _glyph():
    tp = TextPath((0, 0), 'G', size=1., prop=FontProperties(family='DejaVu Sans', weight='bold'))
    polys, cur = [], []
    for curve, code in tp.iter_bezier():              # sample the glyph's Bezier curves finely
        if code == Path.MOVETO:
            if len(cur) > 2:
                polys.append(np.array(cur))
            cur = [curve.control_points[0]]
        else:
            cur.extend(curve(np.linspace(0, 1, 16)[1:]))
    if len(cur) > 2:
        polys.append(np.array(cur))
    v = np.concatenate(polys)
    c, h = (v.min(0) + v.max(0)) / 2, v[:, 1].max() - v[:, 1].min()
    return [Path((p - c) / h) for p in polys]


GLYPH = _glyph()


def in_glyph(n, elev, azim):
    """Whether the unit normals n (world frame) fall inside the G drawn on the front of the sphere."""
    d, right, up = camera(elev, azim)
    a, b, c = n @ right, n @ up, n @ d
    uv = np.c_[np.arctan2(a, c), np.arcsin(np.clip(b, -1, 1))] / G_HEIGHT
    odd = np.zeros(len(n), bool)
    for P in GLYPH:                                   # even-odd rule over the letter's contours
        odd ^= P.contains_points(uv)
    return odd & (c > 0)


def ball_triangles(center, radius, elev, azim, n_lon=240, n_lat=120):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, n_lon + 1), np.linspace(0, np.pi, n_lat + 1))
    S = np.stack([np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v)], -1)
    q = np.stack([S[:-1, :-1], S[:-1, 1:], S[1:, 1:], S[1:, :-1]], 2).reshape(-1, 4, 3)
    unit = np.concatenate([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    n = unit.mean(1)
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    k = .45 + .55 * np.abs(n @ LIGHT)
    rgb = np.where(in_glyph(n, elev, azim)[:, None], WHITE, GREEN) * k[:, None]
    return unit * radius + np.asarray(center), np.c_[np.clip(rgb, 0, 1), np.ones(len(n))]


def ball_image(elev, azim, n=240, ss=3):
    """RGBA image (origin lower) of the ball seen from (elev, azim), supersampled ss x ss."""
    d, right, up = camera(elev, azim)
    m = n * ss
    x = (np.arange(m) + .5) / m * 2 - 1
    X, Y = np.meshgrid(x, x)
    r2 = X ** 2 + Y ** 2
    Z = np.sqrt(np.clip(1 - r2, 0, None))
    N = (X[..., None] * right + Y[..., None] * up + Z[..., None] * d).reshape(-1, 3)
    k = .45 + .55 * np.abs(N @ LIGHT)
    rgb = np.where(in_glyph(N, elev, azim)[:, None], WHITE, GREEN) * k[:, None]
    img = np.c_[np.clip(rgb, 0, 1), (r2 <= 1).reshape(-1)].reshape(m, m, 4)
    img[..., :3] *= img[..., 3:]                      # premultiply, average, un-premultiply
    img = img.reshape(n, ss, n, ss, 4).mean((1, 3))
    img[..., :3] /= np.maximum(img[..., 3:], 1e-9)
    return img
