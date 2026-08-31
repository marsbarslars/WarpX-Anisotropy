"""Find the maximum and minimum coordinates of an STL model."""

from stl import mesh
from pathlib import Path

### USER-DEFINED VARIABLES ###
INPUT_PATH = Path('copyCylinder.stl').resolve()

# Read STL
model = mesh.Mesh.from_file(INPUT_PATH)

# Find bounding box
min_coords = model.vectors.min(axis=(0, 1))
max_coords = model.vectors.max(axis=(0, 1))
size = max_coords - min_coords

print("Bounding box:")
print(f'x_min: {min_coords[0]:.3f}        x_max: {max_coords[0]:.3f}        x_size: {size[0]:.3f}')
print(f'y_min: {min_coords[1]:.3f}        y_max: {max_coords[1]:.3f}        y_size: {size[1]:.3f}')
print(f'z_min: {min_coords[2]:.3f}        z_max: {max_coords[2]:.3f}        z_size: {size[2]:.3f}')