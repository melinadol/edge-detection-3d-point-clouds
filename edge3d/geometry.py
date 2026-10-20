"""Photogrammetric and analytical-geometry building blocks.

Conventions (same as the original thesis code):
- Image coordinates (x, y) are metric (metres) and referred to the principal
  point; the third image coordinate is -c (c = calibrated focal length).
- Exterior orientation: projection centre O = (X0, Y0, Z0) and rotation angles
  omega, phi, kappa. The rotation matrix maps object-space vectors to the
  image frame, so the collinearity equations read

      x = -c * r1 . (P - O) / r3 . (P - O)
      y = -c * r2 . (P - O) / r3 . (P - O)
"""
from __future__ import annotations

import numpy as np

ANGLE_UNITS = {"deg": np.pi / 180.0, "grad": np.pi / 200.0, "rad": 1.0}


def _r_omega(w):
    c, s = np.cos(w), np.sin(w)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])


def _r_phi(f):
    c, s = np.cos(f), np.sin(f)
    return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])


def _r_kappa(k):
    c, s = np.cos(k), np.sin(k)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])


def rotation_matrix(omega, phi, kappa, units="deg", order="kpw"):
    """Rotation matrix R(omega, phi, kappa).

    order="kpw" gives R = R_kappa @ R_phi @ R_omega (the classical
    photogrammetric sequence, used in the thesis simulation).
    order="kwp" gives R = R_kappa @ R_omega @ R_phi, the sequence used in the
    thesis for the Agisoft PhotoScan exterior orientations of the Naxos images
    (angles exported in grads).
    """
    s = ANGLE_UNITS[units]
    w, f, k = omega * s, phi * s, kappa * s
    if order == "kpw":
        return _r_kappa(k) @ _r_phi(f) @ _r_omega(w)
    if order == "kwp":
        return _r_kappa(k) @ _r_omega(w) @ _r_phi(f)
    raise ValueError(f"unknown rotation order {order!r}")


def project(points, c, O, R):
    """Collinearity equations: 3D object points (N,3) -> image coords (N,2)."""
    d = (np.atleast_2d(points) - O) @ R.T  # rows: R @ (P - O)
    return -c * d[:, :2] / d[:, 2:3]


def image_to_object(xy, c, O, R):
    """Rigid-body transformation of image points (x, y, -c) into object space.

    Returns points on the image plane expressed in the object coordinate
    system:  P = O + R^T (x, y, -c).  (Eq. 3 of the paper.)
    """
    xy = np.atleast_2d(xy)
    p = np.column_stack([xy, np.full(len(xy), -c)])
    return O + p @ R  # == (R.T @ p.T).T


def line_2d(p1, p2):
    """Slope/intercept (a, b) of the image line y = a x + b through p1, p2."""
    (x1, y1), (x2, y2) = p1, p2
    a = (y2 - y1) / (x2 - x1)
    return a, y1 - a * x1


def plane_from_points(p1, p2, p3):
    """Plane (A, B, C, D) with Ax + By + Cz + D = 0 through three points."""
    p1, p2, p3 = (np.asarray(p, float).ravel() for p in (p1, p2, p3))
    n = np.cross(p2 - p1, p3 - p1)
    return np.append(n, -n @ p1)


def point_plane_distance(points, plane):
    """Signed perpendicular distance of each point (N,3) from a plane."""
    n = plane[:3]
    return (np.atleast_2d(points) @ n + plane[3]) / np.linalg.norm(n)


def plane_angle(plane1, plane2, units="deg"):
    """Dihedral angle between two planes."""
    n1, n2 = plane1[:3], plane2[:3]
    cos = abs(n1 @ n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
    ang = np.arccos(np.clip(cos, -1.0, 1.0))
    return ang / ANGLE_UNITS[units]
