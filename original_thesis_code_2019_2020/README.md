# Original thesis code (2018–2020)

The code exactly as written during the thesis, kept for the record. Comments and
file names are in Greeklish. For a clean, tested version use the `edge3d` package
in the repository root.

Notes:
- The scripts in `1_simulation/` and `2_error_analysis/*.py` are **Python 2** (Spyder 2);
  the notebooks are Python 3.
- The notebooks expect the Naxos point cloud (`naosDimitrasNaxosMelina_PointcCloud.txt`)
  and hard-coded Windows paths. The data is not included (it belongs to the
  NTUA Laboratory of Photogrammetry), so the notebooks are for reading, not re-running.
- `ransac_etoimos_algorithmos.py` (not included) was the scikit-learn RANSAC example, used for learning.

| Folder / file | What it is |
|---|---|
| `1_simulation/starting.py`, `dokimi*.py` | First experiments (Dec 2018 – Mar 2019): synthetic edge of 30 points, rotation matrix R(ω,φ,κ), collinearity, plane Ax+By+Cz+D=0 through the edge and the projection centre |
| `1_simulation/gewrg_apli_per.py`, `apo_gewrgop_sto_provlima.py` | Simple case and full case of the problem as formulated with the supervisor: image points → object space → plane check |
| `1_simulation/teliko_1o_meros_pseudosystima.py` | Final version of part 1 on the simulated ("pseudo") system |
| `1_simulation/create_random_point_cloud_with_line.ipynb` | Synthetic point cloud (~1000 random points + 200-point edge) used in the paper's simulation |
| `2_error_analysis/epalitheusi_epipedwn_me_alla_simeia_eytheias_KAI_SFALMATA.py` | Verification that all edge points satisfy the plane, and sensitivity of the plane to coordinate errors |
| `2_error_analysis/sfalmata_parametrwn_epipedou_ABCD.ipynb` | Symbolic error propagation (SymPy) for A, B, C, D and the angular error of the plane |
| `2_error_analysis/pollaplasiastes_Lagrange.ipynb` | Lagrange-multiplier formulation experiments |
| `3_naxos_pipeline/project_1_4G0R6724.ipynb` … `project_4_4G0R6615.ipynb` | Full pipeline on the Temple of Demeter (Naxos) for four images: edge plane, iterative RANSAC, acceptance criteria, export to CloudCompare, Dx/Dy statistics |
| `3_naxos_pipeline/Ransac_loops.ipynb`, `ransac_3d_first_attempt.ipynb` | Development of the iterative RANSAC loop |
