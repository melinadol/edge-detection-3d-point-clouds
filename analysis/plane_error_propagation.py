"""How precisely is the edge plane defined by image measurements?

The plane is fixed by the projection centre O and two points of the image edge
transformed into object space. The image points span only a few millimetres,
so the plane could, in principle, be weakly defined. This script propagates
the image-measurement error to the plane parameters (first order + Monte Carlo)
and expresses it as an angular error and as a lateral offset of the plane at
typical object distances.

Run:  python analysis/plane_error_propagation.py
"""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge3d import (image_to_object, plane_from_points,  # noqa: E402
                    plane_orientation_error, plane_parameter_sigma,
                    rotation_matrix)

# Exterior orientation of image 4G0R6724 (Temple of Demeter, Naxos), as used
# in the thesis; angles in grads, PhotoScan rotation order.
O = np.array([1013.2260, 1015.2280, 104.0000])
R = rotation_matrix(-61.4238, 122.3646, 80.5447, units="grad", order="kwp")
c = 0.0245
p1, p2 = np.array([0.0037, 0.0087]), np.array([-0.0091, 0.0085])  # image edge (m)

pixel = 0.036 / 6000                     # 6 µm pixel, full-frame sensor
P1, P2 = image_to_object(np.array([p1, p2]), c, O, R)
plane = plane_from_points(P1, P2, O)

print("Edge plane  A, B, C, D =", np.array2string(plane, precision=6))
for k in (0.5, 1.0, 2.0):
    sigma = k * pixel
    s_abcd, _ = plane_parameter_sigma(P1, P2, O, sigma)
    ang, ang_mc = plane_orientation_error(P1, P2, O, sigma, n_mc=5000)
    print(f"\nimage measurement precision: {k:.1f} px ({1e6 * sigma:.0f} µm)")
    print("  σ(A, B, C, D)       =", np.array2string(s_abcd, precision=8))
    print(f"  angular σ (analytic) = {ang:.2e} rad = {ang * 200 / np.pi * 1e4:.1f} cc")
    print(f"  angular σ (MC)       = {ang_mc:.2e} rad")
    # what matters is how far the perturbed plane passes from the true edge
    rng = np.random.default_rng(1)
    for dist in (5, 10, 20):
        rays = np.array([P1, P2, (P1 + P2) / 2]) - O
        Q = O + rays / np.linalg.norm(rays, axis=1)[:, None] * dist  # edge at `dist`
        off = []
        for _ in range(2000):
            q1, q2 = image_to_object(np.array([p1, p2]) + rng.normal(0, sigma, (2, 2)),
                                     c, O, R)
            pl = plane_from_points(q1, q2, O)
            off.append((Q @ pl[:3] + pl[3]) / np.linalg.norm(pl[:3]))
        rms = 1000 * np.sqrt(np.mean(np.square(off)))
        print(f"  RMS offset of the plane from an edge at {dist:2d} m: {rms:5.2f} mm "
              f"(≈ σ·d/c = {1000 * sigma * dist / c:4.2f} mm)")
