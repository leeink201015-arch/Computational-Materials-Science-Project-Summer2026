import numpy as np
from scipy.spatial.transform import Rotation as R

tilt_angles = np.array([0, 30, 60, 90, 120, 150, 180])
azimuth_angles = np.array([0, 60, 120, 180, 240, 300])

atom_coordinate_original = np.array([
    [4.0, 4.0, 2.0],
    [4.0, 4.0, 4.0],
    [4.0, 4.0, 6.0]
])

rotation_centre = atom_coordinate_original[1]
relative_coordinates=(
    atom_coordinate_original - rotation_centre
)

# We need to make sure that we are rotating around the point [4,4,4] not the origin
for tilt_angle in tilt_angles:
    R_tilted = np.array([
        [1, 0, 0],
        [0, np.cos(np.radians(tilt_angle)), -np.sin(np.radians(tilt_angle))],
        [0, np.sin(np.radians(tilt_angle)), np.cos(np.radians(tilt_angle))]
    ])

    for azimuth_angle in azimuth_angles:
        R_azimuth = np.array([
            [np.cos(np.radians(azimuth_angle)), -np.sin(np.radians(azimuth_angle)), 0],
            [np.sin(np.radians(azimuth_angle)), np.cos(np.radians(azimuth_angle)), 0],
            [0, 0, 1]
        ])

        # Matrix multiplication (it's not communtative), Final Point P_f = R_azimuth x (R_tilt x P_i) (I'll convert this to LaTeX later whatever)

        atom_coordinate_tilted = atom_coordinate_original @ R_tilted.T
        atom_coordinate_rotated = atom_coordinate_tilted @ R_azimuth.T

        atom_coordinate_final=(
            atom_coordinate_rotated + rotation_centre
        )

        print(atom_coordinate_final)