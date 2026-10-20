"""edge3d – 3D edge detection in point clouds using oriented digital images.

Refactored implementation of the NTUA Dipl.-Eng. thesis of M.M. Dolapsaki (2020)
and of Dolapsaki & Georgopoulos, ISPRS Int. J. Geo-Inf. 2021, 10, 229.
"""
from .geometry import (rotation_matrix, project, image_to_object, line_2d,
                       plane_from_points, point_plane_distance, plane_angle)
from .pipeline import Camera, detect_edge, ransac_line, save_xyz
from .errors import plane_parameter_sigma, plane_orientation_error

__all__ = [
    "rotation_matrix", "project", "image_to_object", "line_2d",
    "plane_from_points", "point_plane_distance", "plane_angle",
    "Camera", "detect_edge", "ransac_line", "save_xyz",
    "plane_parameter_sigma", "plane_orientation_error",
]
__version__ = "1.0.0"
