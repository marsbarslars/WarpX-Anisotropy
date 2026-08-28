from stl import mesh

# Input/output files
input_file = "w7x.stl"
output_file = "w7x_scaled.stl"

# Scale factor
scale = 0.5

# Read STL
stl_mesh = mesh.Mesh.from_file(input_file)

# Scale all vertex coordinates
stl_mesh.vectors *= scale

# Save scaled STL
stl_mesh.save(output_file)

print(f"Saved scaled STL to {output_file}")