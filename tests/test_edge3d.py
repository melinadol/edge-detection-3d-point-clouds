from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge3d import (Camera, detect_edge, image_to_object, plane_from_points,  # noqa: E402
                    point_plane_distance, project, rotation_matrix)


def test_rotation_is_orthonormal():
    for order in ("kpw", "kwp"):
        R = rotation_matrix(12.3, -45.6, 78.9, units="grad", order=order)
        assert np.allclose(R @ R.T, np.eye(3))
        assert np.isclose(np.linalg.det(R), 1.0)


def test_collinearity_roundtrip():
    R = rotation_matrix(88, 3, -2)
    O = np.array([1.0, -9.0, 1.5])
    c = 0.0245
    P = np.array([[0.5, 1.0, 2.0], [3.0, 2.0, 0.2]])
    xy = project(P, c, O, R)
    E = image_to_object(xy, c, O, R)       # points on the image plane in object space
    for p, e in zip(P, E):                 # O, e, p must be collinear
        assert np.linalg.norm(np.cross(e - O, p - O)) < 1e-9


def test_edge_points_lie_on_edge_plane():
    R = rotation_matrix(88, 3, -2)
    O = np.array([4.0, -9.0, 1.5])
    c = 0.0245
    edge = np.array([[1.0, 0.0, 0.5], [6.0, 0.0, 2.5]])
    p1, p2 = project(edge, c, O, R)
    E1, E2 = image_to_object(np.array([p1, p2]), c, O, R)
    plane = plane_from_points(E1, E2, O)
    pts = edge[0] + np.linspace(0, 1, 20)[:, None] * (edge[1] - edge[0])
    assert np.abs(point_plane_distance(pts, plane)).max() < 1e-9


def test_detect_edge_on_simple_scene():
    rng = np.random.default_rng(0)
    # a vertical wall (y = 0) with a box in front; detect the box's front-top edge
    wall = np.column_stack([rng.uniform(0, 6, 40000), np.zeros(40000),
                            rng.uniform(0, 3, 40000)])
    u = np.arange(1, 5, 0.01)
    v = np.arange(0, 1, 0.01)
    U, V = np.meshgrid(u, v)
    front = np.column_stack([U.ravel(), np.full(U.size, -1.0), V.ravel()])
    W, Y = np.meshgrid(u, np.arange(-1, 0, 0.01))
    top = np.column_stack([W.ravel(), Y.ravel(), np.ones(W.size)])
    cloud = np.vstack([wall, front, top]) + rng.normal(0, 0.001, (len(wall) + len(front) + len(top), 3))

    cam = Camera(c=0.0245, O=[3.0, -10.0, 2.5], R=rotation_matrix(80, 0, 0))
    a, b = np.array([1.5, -1.0, 1.0]), np.array([4.5, -1.0, 1.0])
    p1, p2 = project(np.array([a, b]), cam.c, cam.O, cam.R)
    res = detect_edge(cloud, p1, p2, cam, image_threshold=2 * cam.frame[0] / 6000)
    assert res.edge is not None
    d = res.edge - a
    dirn = (b - a) / np.linalg.norm(b - a)
    off = np.linalg.norm(d - np.outer(d @ dirn, dirn), axis=1)
    assert np.sqrt(np.mean(off ** 2)) < 0.01
