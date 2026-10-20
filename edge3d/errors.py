"""Error propagation for the edge plane (thesis chapter on accuracy evaluation).

The plane is defined by three points (two edge points transformed from the
image into object space and the projection centre). Assuming independent
coordinate errors with standard deviation `sigma`, the first-order
propagation  Σ_ABCD = σ² J Jᵀ  gives the precision of the plane parameters,
and the angular precision of the plane follows from the unit normal.
"""
from __future__ import annotations

import numpy as np

from .geometry import plane_from_points


def _jacobian(f, x, h=1e-7):
    x = np.asarray(x, float)
    f0 = f(x)
    J = np.empty((f0.size, x.size))
    for i in range(x.size):
        dx = np.zeros_like(x)
        dx[i] = h * max(1.0, abs(x[i]))
        J[:, i] = (f(x + dx) - f(x - dx)) / (2 * dx[i])
    return J


def plane_parameter_sigma(p1, p2, p3, sigma):
    """Standard deviations of (A, B, C, D) and their covariance matrix."""
    x = np.concatenate([np.ravel(p1), np.ravel(p2), np.ravel(p3)])
    J = _jacobian(lambda v: plane_from_points(v[:3], v[3:6], v[6:]), x)
    cov = sigma ** 2 * J @ J.T
    return np.sqrt(np.diag(cov)), cov


def plane_orientation_error(p1, p2, p3, sigma, n_mc=0, rng=0):
    """Angular standard error (rad) of the plane normal.

    Analytical first-order value; if n_mc > 0 a Monte-Carlo estimate is
    returned as well (RMS angle between perturbed and true planes).
    """
    x = np.concatenate([np.ravel(p1), np.ravel(p2), np.ravel(p3)])

    def unit_normal(v):
        n = plane_from_points(v[:3], v[3:6], v[6:])[:3]
        return n / np.linalg.norm(n)

    J = _jacobian(unit_normal, x)
    cov_n = sigma ** 2 * J @ J.T
    sigma_ang = float(np.sqrt(np.trace(cov_n)))
    if not n_mc:
        return sigma_ang
    rng = np.random.default_rng(rng)
    n0 = unit_normal(x)
    ang = []
    for _ in range(n_mc):
        n = unit_normal(x + rng.normal(0, sigma, x.size))
        ang.append(np.arccos(np.clip(abs(n @ n0), -1, 1)))
    return sigma_ang, float(np.sqrt(np.mean(np.square(ang))))
