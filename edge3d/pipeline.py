"""3D edge detection in a point cloud from an edge measured on an oriented image.

Part 1 – plane:  the edge measured on the image and the projection centre
                 define a plane in object space. Every point of the 3D edge
                 lies on it, so the points of the cloud close to that plane
                 are kept.
Part 2 – edge:   the kept points contain several linear concentrations (the
                 plane cuts many surfaces). Lines are extracted iteratively
                 with RANSAC and each candidate is tested with acceptance
                 criteria:
                   1. back-projected onto the image it must fall inside the
                      physical image frame and inside the extent of the
                      measured 2D edge (parts outside are discarded);
                   2. back-projected, it must agree with the image line
                      (|dx|, |dy| below a threshold);
                   3. among the accepted candidates the one nearest to the
                      projection centre is selected (it is the visible one;
                      edges further away along the same plane are occluded).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .geometry import (image_to_object, line_2d, plane_from_points,
                       point_plane_distance, project)


@dataclass
class Camera:
    c: float                      # calibrated focal length (m)
    O: np.ndarray                 # projection centre (X0, Y0, Z0)
    R: np.ndarray                 # 3x3 rotation matrix
    frame: tuple = (0.036, 0.024)  # sensor size (m) – full-frame DSLR

    def __post_init__(self):
        self.O = np.asarray(self.O, float)


@dataclass
class Candidate:
    points: np.ndarray
    dx: np.ndarray
    dy: np.ndarray
    perp: np.ndarray              # perpendicular image distance from the 2D edge
    distance: float               # mean distance from the projection centre
    accepted: bool


@dataclass
class EdgeResult:
    plane: np.ndarray
    plane_points: np.ndarray
    candidates: list = field(default_factory=list)
    edge: np.ndarray | None = None


def ransac_line(points, threshold, max_trials=1000, min_inliers=10, rng=None):
    """Fit a 3D line with RANSAC. Returns (point, direction, inlier_mask)."""
    rng = np.random.default_rng(rng)
    n = len(points)
    best = None
    best_count = 0
    for _ in range(max_trials):
        i, j = rng.choice(n, 2, replace=False)
        d = points[j] - points[i]
        norm = np.linalg.norm(d)
        if norm == 0:
            continue
        d /= norm
        v = points - points[i]
        resid = np.linalg.norm(v - np.outer(v @ d, d), axis=1)
        mask = resid < threshold
        cnt = mask.sum()
        if cnt > best_count:
            best_count, best = cnt, mask
    if best is None or best_count < min_inliers:
        return None
    # least-squares refinement on the inliers (principal axis)
    inl = points[best]
    centre = inl.mean(axis=0)
    direction = np.linalg.svd(inl - centre)[2][0]
    v = points - centre
    resid = np.linalg.norm(v - np.outer(v @ direction, direction), axis=1)
    return centre, direction, resid < threshold


def inside_frame(points, cam: Camera):
    xy = project(points, cam.c, cam.O, cam.R)
    hw, hh = cam.frame[0] / 2, cam.frame[1] / 2
    return (np.abs(xy[:, 0]) <= hw) & (np.abs(xy[:, 1]) <= hh)


def inside_segment(points, cam: Camera, p1, p2, margin=0.02):
    """Points whose image falls within the measured segment p1–p2 (± margin)."""
    xy = project(points, cam.c, cam.O, cam.R)
    p1, p2 = np.asarray(p1, float), np.asarray(p2, float)
    d = p2 - p1
    t = (xy - p1) @ d / (d @ d)
    return (t >= -margin) & (t <= 1 + margin)


def backprojection_residuals(points, cam: Camera, a, b):
    """dx, dy of back-projected points from the image line y = a x + b (Eq. 4–5)."""
    xy = project(points, cam.c, cam.O, cam.R)
    x, y = xy[:, 0], xy[:, 1]
    dy = np.abs(y - (a * x + b))
    dx = np.abs(x - (y - b) / a) if a != 0 else np.full_like(x, np.inf)
    return dx, dy


def detect_edge(cloud, p1, p2, cam: Camera, plane_threshold=0.006,
                line_threshold=0.006, image_threshold=None, max_lines=10,
                min_inliers=30, rng=0):
    """Detect the 3D edge of `cloud` that is seen on the image as p1–p2.

    cloud            (N,3) object coordinates
    p1, p2           two image points (m) on the measured 2D edge
    plane_threshold  max distance of a cloud point from the edge plane (m)
    line_threshold   RANSAC residual threshold for 3D lines (m)
    image_threshold  max |dx|,|dy| on the image (m); default 1 pixel of a
                     36 mm / 6000 px sensor = 6 µm
    """
    if image_threshold is None:
        image_threshold = cam.frame[0] / 6000
    a, b = line_2d(p1, p2)

    # ---- Part 1: plane through the projection centre and the image edge
    P1, P2 = image_to_object(np.array([p1, p2]), cam.c, cam.O, cam.R)
    plane = plane_from_points(P1, P2, cam.O)
    dist = np.abs(point_plane_distance(cloud, plane))
    on_plane = cloud[dist <= plane_threshold]
    # keep only points in front of the camera
    depth = (on_plane - cam.O) @ cam.R.T
    on_plane = on_plane[depth[:, 2] < 0]

    result = EdgeResult(plane=plane, plane_points=on_plane)

    # ---- Part 2: iterative RANSAC + acceptance criteria
    remaining = on_plane.copy()
    for _ in range(max_lines):
        if len(remaining) < min_inliers:
            break
        fit = ransac_line(remaining, line_threshold, min_inliers=min_inliers,
                          rng=rng)
        if fit is None:
            break
        _, _, mask = fit
        pts = remaining[mask]
        remaining = remaining[~mask]

        pts = pts[inside_frame(pts, cam)                # criterion 1
                  & inside_segment(pts, cam, p1, p2)]
        if len(pts) < min_inliers:
            continue
        dx, dy = backprojection_residuals(pts, cam, a, b)
        # criterion 2 – perpendicular image distance (= dy·cos(atan a)),
        # robust to the noise of the scanner points via the median
        perp = dy / np.sqrt(1 + a * a)
        accepted = bool(np.median(perp) <= image_threshold)
        d = float(np.linalg.norm(pts - cam.O, axis=1).mean())
        result.candidates.append(Candidate(pts, dx, dy, perp, d, accepted))

    accepted = [cnd for cnd in result.candidates if cnd.accepted]
    if accepted:                                          # criterion 3
        result.edge = min(accepted, key=lambda cnd: cnd.distance).points
    return result


def save_xyz(points, path):
    """Export points as an ASCII XYZ file (opens directly in CloudCompare)."""
    np.savetxt(path, points, fmt="%.4f")
