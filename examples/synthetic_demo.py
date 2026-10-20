"""End-to-end demo on a synthetic monument-like scene.

A stepped block (two steps + wall) and a column standing behind it are sampled as a noisy "laser-scanner" point
cloud. A virtual full-frame camera with known interior/exterior orientation
looks at it. One straight edge is measured on the image (with half a pixel of
measurement noise) and the algorithm recovers the corresponding 3D points.

Run:  python examples/synthetic_demo.py
"""
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge3d import Camera, detect_edge, project, rotation_matrix  # noqa: E402

FIG = Path(__file__).resolve().parents[1] / "docs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)


def sample_box(x0, x1, y0, y1, z0, z1, spacing, noise):
    """Points on the 5 visible faces of an axis-aligned box (no bottom)."""
    pts = []

    def grid(u0, u1, v0, v1):
        u = np.arange(u0, u1, spacing)
        v = np.arange(v0, v1, spacing)
        return np.meshgrid(u, v)

    for y in (y0, y1):                         # front / back
        X, Z = grid(x0, x1, z0, z1)
        pts.append(np.column_stack([X.ravel(), np.full(X.size, y), Z.ravel()]))
    for x in (x0, x1):                         # sides
        Y, Z = grid(y0, y1, z0, z1)
        pts.append(np.column_stack([np.full(Y.size, x), Y.ravel(), Z.ravel()]))
    X, Y = grid(x0, x1, y0, y1)                # top
    pts.append(np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, z1)]))
    P = np.vstack(pts)
    return P + rng.normal(0, noise, P.shape)


def main():
    spacing, noise = 0.02, 0.002                       # 2 cm spacing, 2 mm noise
    cloud = np.vstack([
        sample_box(0.0, 8.0, 0.0, 3.0, 0.0, 0.4, spacing, noise),   # step 1
        sample_box(0.3, 7.7, 0.3, 2.7, 0.4, 0.8, spacing, noise),   # step 2
        sample_box(0.8, 7.2, 0.8, 2.2, 0.8, 4.5, spacing, noise),   # wall
        # a column behind the monument, hidden from the camera; one of its
        # vertical edges lies (almost) on the same plane as the measured edge
        sample_box(-0.66, -0.16, 5.0, 5.5, 0.0, 6.0, spacing, noise),
    ])
    print(f"point cloud: {len(cloud):,} points")

    # camera: full-frame DSLR, c = 24.5 mm, ~9 m in front of the monument
    cam = Camera(c=0.0245, O=np.array([4.2, -9.0, 1.6]),
                 R=rotation_matrix(88.0, 3.0, -2.0, units="deg"))

    # true edge: left vertical edge of the wall (x = 0.8, y = 0.8)
    true_a, true_b = np.array([0.8, 0.8, 1.0]), np.array([0.8, 0.8, 4.3])
    p1, p2 = project(np.array([true_a, true_b]), cam.c, cam.O, cam.R)
    px = cam.frame[0] / 6000
    p1 = p1 + rng.normal(0, 0.5 * px, 2)                 # measurement noise
    p2 = p2 + rng.normal(0, 0.5 * px, 2)

    res = detect_edge(cloud, p1, p2, cam, plane_threshold=0.006,
                      line_threshold=0.006, image_threshold=1.0 * px)

    print(f"points within 6 mm of the edge plane: {len(res.plane_points):,}")
    for i, cnd in enumerate(res.candidates, 1):
        med = np.median(cnd.perp) / px
        print(f"  candidate #{i}: {len(cnd.points):4d} pts, "
              f"mean dist to O = {cnd.distance:6.2f} m, "
              f"median image offset = {med:7.2f} px, accepted = {cnd.accepted}")

    edge = res.edge
    d = edge - true_a
    u = (true_b - true_a) / np.linalg.norm(true_b - true_a)
    off = np.linalg.norm(d - np.outer(d @ u, u), axis=1)
    print(f"\nselected edge: {len(edge)} points, distance from the true edge: "
          f"RMS {1000 * np.sqrt(np.mean(off**2)):.1f} mm, "
          f"max {1000 * off.max():.1f} mm")

    # ---------- figure 1: 3D view with the edge plane
    from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    def status(cnd):
        if cnd.points is edge:
            return "selected"
        return "accepted, occluded (further)" if cnd.accepted else "rejected"

    style = {"selected": ("crimson", 7), "accepted, occluded (further)": ("#2b6cb0", 3),
             "rejected": ("#9aa5b1", 2)}
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    sub = cloud[rng.choice(len(cloud), 30000, replace=False)]
    ax.scatter(*sub.T, s=0.3, c="0.72", alpha=0.35)
    for cnd in res.candidates:
        col, s = style[status(cnd)]
        ax.scatter(*cnd.points.T, s=s, color=col)
    # plane wedge: projection centre -> rays through the image edge end points
    from edge3d import image_to_object
    E = image_to_object(np.array([p1, p2]), cam.c, cam.O, cam.R)
    far = cam.O + (E - cam.O) * (16.5 / np.linalg.norm(E - cam.O, axis=1))[:, None]
    ax.add_collection3d(Poly3DCollection([[cam.O, far[0], far[1]]],
                                         alpha=0.12, color="crimson"))
    for f in far:
        ax.plot(*np.array([cam.O, f]).T, color="crimson", lw=0.6, alpha=0.6)
    ax.scatter(*cam.O, marker="^", s=40, c="k")
    ax.set_xlabel("X (m)"); ax.set_ylabel("Y (m)"); ax.set_zlabel("Z (m)")
    ax.view_init(elev=24, azim=-58)
    ax.set_box_aspect((9, 15, 6.5))
    handles = [Line2D([], [], ls="", marker="o", color=style[k][0], label=k)
               for k in style] + [
        Line2D([], [], ls="", marker="^", color="k", label="projection centre O"),
        Line2D([], [], color="crimson", alpha=0.5, lw=6, label="edge plane (O + image edge)")]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=False)
    ax.set_title("The plane through O and the image edge cuts the cloud;\n"
                 "candidate 3D lines are tested and the nearest valid one is kept",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "synthetic_3d.png", dpi=150)

    # ---------- figure 2: acceptance criteria per candidate
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for i, cnd in enumerate(res.candidates, 1):
        col, _ = style[status(cnd)]
        y = np.median(cnd.perp) / px
        ax.scatter(cnd.distance, y, s=20 + len(cnd.points) / 2, color=col,
                   edgecolor="white", zorder=3)
        ax.annotate(f"#{i}", (cnd.distance, y), xytext=(8, 4),
                    textcoords="offset points", fontsize=9)
    ax.axhline(1.0, color="0.3", ls="--", lw=1)
    ax.text(ax.get_xlim()[0], 1.03, " 1-pixel threshold", va="bottom", fontsize=8,
            color="0.3")
    ax.set_xlabel("mean distance of candidate from projection centre (m)")
    ax.set_ylabel("median image offset (pixels)")
    ax.legend(handles=handles[:3], fontsize=8, frameon=False, loc="upper right")
    ax.set_title("Candidates: below the threshold = accepted; nearest accepted = edge",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "synthetic_candidates.png", dpi=150)
    print(f"figures written to {FIG}")


if __name__ == "__main__":
    main()
