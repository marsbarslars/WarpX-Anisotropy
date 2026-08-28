import h5py
import numpy as np

filename = "example-femm-3d.h5"

with h5py.File(filename, "r+") as f:

    # # Read
    # x = f["x"][:]
    # y = f["y"][:]
    # z = f["z"][:]

    # Look at the B mesh
    data = f['data']
    print('Attributes of data:')
    for name, value in data.attrs.items():
        print(f"{name}: {value}")

    mesh = f['data/1/meshes']
    print('\nAttributes of mesh:')
    for name, value in mesh.attrs.items():
        print(f"{name}: {value}")

    B = f["data/1/meshes/B"]

    print("\nAttributes of B:")
    for name, value in B.attrs.items():
        print(f"{name}: {value}")

    print("\nAttributes of B/x:")
    for name, value in B["x"].attrs.items():
        print(f"{name}: {value}")
    # By = f["By"][:]
    # Bz = f["Bz"][:]

    # # Alter
    # Bx *= 1.1
    # By *= 1.1
    # Bz *= 1.1

    # # Write altered values back
    # f["Bx"][:] = Bx
    # f["By"][:] = By
    # f["Bz"][:] = Bz

    # # Add a new dataset
    # Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)

    # f.create_dataset("Bmag", data=Bmag)

    print("Done.")

    def print_structure(name, obj):
        print(name)

    f.visititems(print_structure)

    dx, dy, dz = f["data/1/meshes/B"].attrs["gridSpacing"]
    x0, y0, z0 = f["data/1/meshes/B"].attrs["gridGlobalOffset"]
    Bx = f['data/1/meshes/B/x'][:]

    Nx, Ny, Nz = Bx.shape

    x = x0 + np.arange(Nx) * dx
    y = y0 + np.arange(Ny) * dy
    z = z0 + np.arange(Nz) * dz

print("Grid shape:", Bx.shape)
print("x:", x)
print("y:", y)
print("z:", z)