"""
Modifies an STL file. Order of steps:
Scale
Mirror
Rotate
Translate
"""

from stl import mesh
from pathlib import Path
import numpy as np

### USER-DEFINED VARIABLES ###
INPUT_PATH = Path('cylinder.stl').resolve()
OUTPUT_PATH = Path('copyCylinder.stl').resolve()
SCALE = 0.001
MIRROR = False
MIRROR_PLANE = 'XY'  # Options: 'XY', 'XZ', 'YZ'
ROTATION_POINT = np.array([0., 0., 0.])
ROTATION_AXIS = np.array([0., 0., 1.])
ROTATION_ANGLE = 0.  # In degrees
TRANSLATION = np.array([0., 0., 0.])


def main() -> None:

    # Throw errors
    if SCALE <= 0:
        raise ValueError("SCALE must be greater than zero")
    if INPUT_PATH == OUTPUT_PATH:
        raise ValueError("input and output paths must be different")
    if OUTPUT_PATH.exists():
        raise FileExistsError(f"output file already exists: {OUTPUT_PATH}")
    if MIRROR_PLANE not in ['XY', 'XZ', 'YZ']:
        raise ValueError("MIRROR_PLANE must be 'XY', 'XZ', or 'YZ'")
    if np.linalg.norm(ROTATION_AXIS) == 0:
        raise ValueError("ROTATION_AXIS must be a non-zero vector")

    # Read STL
    model = mesh.Mesh.from_file(INPUT_PATH)

    # Scale all vertex coordinates
    model.vectors *= SCALE

    # Mirror the model if requested
    if MIRROR == True:
        model.vectors = model.vectors[:, ::-1, :]
        if MIRROR_PLANE == 'XY':
            model.vectors[:, :, 2] *= -1
        elif MIRROR_PLANE == 'XZ':
            model.vectors[:, :, 1] *= -1
        elif MIRROR_PLANE == 'YZ':
            model.vectors[:, :, 0] *= -1

    # Rotate the model about an axis coming from a specific point
    axis = ROTATION_AXIS / np.linalg.norm(ROTATION_AXIS)
    theta = np.radians(ROTATION_ANGLE)
    K = np.array([
        [0,        -axis[2],  axis[1]],
        [axis[2],   0,       -axis[0]],
        [-axis[1],  axis[0],  0]
    ])
    R = (
        np.eye(3)
        + np.sin(theta) * K
        + (1 - np.cos(theta)) * (K @ K)
    )
    model.vectors = (
        (model.vectors - ROTATION_POINT) @ R.T
        + ROTATION_POINT
    )

    # Translate the model
    model.vectors += TRANSLATION

    # Save scaled STL
    model.save(OUTPUT_PATH)

    print(f"Created: {OUTPUT_PATH}")
    print(f"Scale Factor: {SCALE}")
    print(f"Mirrored: {MIRROR}")
    if MIRROR == True:
        print(f"Mirror Plane: {MIRROR_PLANE}")
    print(f"Rotation Point: {ROTATION_POINT}")
    print(f"Rotation Axis: {ROTATION_AXIS}")
    print(f"Rotation Angle: {ROTATION_ANGLE} degrees")
    print(f"Translation: {TRANSLATION}")

if __name__ == "__main__":
    main()