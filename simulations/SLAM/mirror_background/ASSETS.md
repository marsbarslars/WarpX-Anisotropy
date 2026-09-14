# Magnetic field and vessel assets

Audited locally on 2026-09-10 using read-only HDF5 and STL inspection. The new workflow references existing assets by path; it does not copy, translate, rescale, or modify the originals.

## Available magnetic field

The default field is the existing tutorial map:

```text
Windows:
C:\Users\ryanv\Desktop\BOOGERAIDS\WIPPL\hackathon-beam-anisotropy\warpx\official_magnetic_mirror\example-femm-3d.h5

WSL:
/mnt/c/Users/ryanv/Desktop/BOOGERAIDS/WIPPL/hackathon-beam-anisotropy/warpx/official_magnetic_mirror/example-femm-3d.h5

SHA256:
f6845dd8046eed415face768abe98748bb624523ea4165ffe0f026d2498245eb

Size: 3,561,320 bytes
```

This is an openPMD 1.1.0 HDF5 file with a Cartesian field mesh at iteration 1. Its metadata identifies openPMD-api 0.15.0 and a creation date of 2023-05-23. The magnetic components are datasets under `/data/1/meshes/B/{x,y,z}`. Each component uses SI units (tesla), and all sampled values are finite. The file also contains constant-zero electric components; the workflow reads the prescribed magnetic field and separately configures its electrostatic model.

| Property | Audited value |
| --- | --- |
| Axis labels and array dimensions | x, y, z; 47 by 47 by 47 |
| Grid spacing | 0.05, 0.05, 0.125 m |
| First node | (-1.15, -1.15, -0.375) m |
| Last node | (1.15, 1.15, 5.375) m |
| Full sampled map minimum magnitude | 0.000124238142 T |
| Full sampled map maximum magnitude | 0.009238652675 T |
| Maximum sampled magnitude inside r < 0.9 m, 0 <= z <= 5 m | 0.006467056685 T |
| On-axis minimum in 0 <= z <= 5 m | 0.001487675040 T |
| On-axis maximum in 0 <= z <= 5 m | 0.004625954116 T |
| On-axis mirror ratio, maximum/minimum | 3.109519211 |

The field covers the whole configured simulation box, x/y from -1 to 1 m and z from 0 to 5 m. The magnitude maxima above are maxima of the stored samples, not estimates from a coil model. The full-map maximum supplies a conservative magnetic-field value for the electron timestep check.

At the center of Lars's injection window, (-0.2, 0, 1) m, trilinear interpolation gives B approximately (-0.000238203, 0, 0.003381358) T and magnitude 0.003389738 T. This window lies within both the field mesh and the vessel bore.

The on-axis ratio of approximately 3.1095 describes this tutorial field. It is not a universal loss-cone threshold for every off-axis birth point or a verified SLAM equilibrium parameter.

## Identity and other copies

The same SHA256 was measured for the copies in the legacy `warpx/mirror_collisions`, `warpx/mirror_angle_scan`, `warpx/mirror_3d_run`, and `Plasma-Hackathon/runs/magnetic-mirror` directories, and for:

```text
/home/ryanv/src/openPMD-example-datasets/example-femm-3d.h5
```

The similarly named file under `Plasma-Hackathon/runs/mirror-validation` is a flattened symbolic-link text file containing `../magnetic-mirror/example-femm-3d.h5`; it is not itself an HDF5 dataset and should not be selected as a field input.

No HDF5 file was found in the shared `WarpX-Anisotropy` checkout, including ignored files. Lars's checked-in Mirror deck requests the filename `example-femm-3d.h5`, which matches the available tutorial asset. We have not established that this hash is the exact copy Lars used in his demonstrated run. The current package records the selected map as a known local candidate rather than claiming that identity.

The shared `simulations/Mirror/hdf5_modify.py` currently specifies a -2.5 m z translation into a separate `copy.h5`. This workflow does not execute that script. Such a translated field would have different coverage and would require a corresponding geometry review; the original map already covers the current z = 0 to 5 m domain.

## Periodic seam limitation

The original field is not an ordinary periodic vector field across z = 0 and z = 5 m. Its transverse components reverse approximately across the ends. Within the r < 0.9 m bore, the measured maximum absolute endpoint differences are:

| Component | Maximum absolute endpoint difference |
| --- | --- |
| Bx | 0.00267927148 T |
| By | 0.00267927148 T |
| Bz | 0.0000182814591 T |

Ordinary periodic particle wrapping can therefore bring a particle into a different local magnetic vector immediately after crossing. Short startup runs with particles initially away from the seam do not validate long-time circulation or mirror end losses. The intended custom boundary and its unresolved provenance are documented in `BOUNDARY_STATUS.md` beside this file. The map is left unchanged.

## Original STL vessel

```text
C:\Users\ryanv\Desktop\BOOGERAIDS\WIPPL\WarpX-Anisotropy\simulations\Mirror\cylinder_OGMirror.stl

SHA256:
75598bc031626a734f58e83b0c7505a3015b9f8afd1e49ad15cc686ceb5479f7
```

The STL is an annular solid shell: inner radius 0.9 m, outer radius 1.0 m, and z extent 0 to 5 m. The annular end faces close the solid wall while leaving the central bore open. Its sampled bounds are x = [-1, 1] m, y = [-0.99909896, 0.99909896] m, and z = [0, 5] m. The slight y-bound difference follows from polygonal sampling of the circle.

Inspection found 1,152 triangles, 576 unique vertices, and 1,728 unique edges, each shared by two triangles. Normals point into the bore on the inner wall and outward at the outer wall, consistent with the solid shell orientation. Its signed enclosed solid volume is approximately 2.982725 m3. These mesh checks do not establish that every downstream embedded-boundary solver setup will converge.

The neighboring `cylinder.stl` has SHA256 `231b46290b00867e97ccbc976d00de332b414adcdb17e5bab7161e6fd73317610`. Its coordinates are 200 times larger: inner/outer radii 180/200 and length 1,000 in raw STL units. The checked-in `stl_modify.py` applies a scale of 0.005 to produce the OG-mirror dimensions. The workflow does not execute that conversion or substitute the unscaled file.

## Wall used by the new background runs

During integration of this package, the electrostatic Poisson solve stalled with the original STL wall, while a grounded analytic bore of radius 0.9 m passed the strict solve used for the startup checks. The cause of the STL solve behavior has not been isolated; this is not a claim that the triangulation is defective.

The package therefore defaults to `--wall analytic-bore`, with:

```text
warpx.eb_implicit_function = "x*x+y*y-0.9*0.9"
warpx.eb_potential(x,y,z,t) = "0.0"
```

This represents the original shell's inner cylindrical surface by a smooth analytic embedded boundary, grounded at 0 V, with an absorbing particle wall. It preserves the intended 0.9 m bore but is an explicit geometry approximation to the STL case. The original STL remains available through `--wall stl` for targeted troubleshooting; it is not the geometry used for the successful default startup runs.
