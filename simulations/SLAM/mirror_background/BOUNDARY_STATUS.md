# Mirror boundary and solver status

Inspected on 2026-09-10. This records local source and public repository evidence; it does not claim that Lars's custom boundary has been reproduced.

## What the shared input requests

`simulations/Mirror/inputs_mirror.txt` sets:

```text
boundary.field_hi = pec pec periodic
boundary.field_lo = pec pec periodic
boundary.particle_hi = absorbing absorbing Periodic_ReflectParVel
boundary.particle_lo = absorbing absorbing Periodic_ReflectParVel
boundary.particle_eb = absorbing
```

Lars introduced these particle boundary settings in commit `b8acce9b91d65b0ddc5102550e8e29c43ffacc02` on 2026-09-04. The source implementing `Periodic_ReflectParVel` is not included in the shared repository.

The current deck also sets `warpx.do_electrostatic = labframe` and `protons.do_not_deposit = 0`. As written, the protons deposit charge and an electrostatic solver computes their electric field, in addition to the prescribed external magnetic field. `protons.initialize_self_fields = 0` does not disable subsequent deposition or the electrostatic solve. There is no explicit Coulomb collision block in the original deck. These are statements about the checked-in input, not proof of which input produced a particular video from Lars.

## Installed source and build

The installed WSL source is `/home/ryanv/src/warpx`, a clean checkout of upstream `https://github.com/BLAST-WarpX/warpx.git`, branch `development`, commit `bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc`. Its only listed local branch is `development`; the remote-tracking branches are `origin/development` and `origin/HEAD`.

The cache at `/home/ryanv/src/warpx/build_gpu_py/CMakeCache.txt` records `WarpX_APP=ON`, `WarpX_COMPUTE=CUDA`, `WarpX_DIMS=3`, and that same source directory. The executable is `/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d`. Generated build metadata records WarpX `bd8c12ff8fc8` and AMReX `26.08-18-g057940244648`; the AMReX checkout is clean at `057940244648b82908cfc486f07ec796bba2b07f`.

Relevant source evidence:

- `Source/Utils/WarpXAlgorithmSelection.H`, lines 144-151: the particle boundary enum contains `Absorbing`, `Open`, `Reflecting`, `Periodic`, `Thermal`, and `None`; it does not contain `Periodic_ReflectParVel`.
- `Source/Particles/ParticleBoundaries.cpp`, lines 77-85: explicitly supplied boundaries are parsed as those enum values. Lines 38-46 require both particle boundaries to be periodic when the corresponding field boundary is periodic.
- `Source/Particles/WarpXParticleContainer.cpp`, lines 2945-2946: periodic wrapping is handled by AMReX.
- `build_gpu_py/_deps/fetchedamrex-src/Src/Particle/AMReX_ParticleUtil.H`, lines 423-457: `enforcePeriodic` adds or subtracts a domain length from position. It does not change momentum.
- `Source/Particles/ParticleBoundaries_K.H`, lines 42-45 and 65-68: ordinary reflecting boundaries reflect the position at the same wall and flag a change of sign of normal momentum.

Therefore the original custom boundary value is unsupported by this installed source. An original-deck run has not been launched for this audit; rejection of the unrecognized enum is a source-level conclusion. A case using standard `periodic` is a different, explicitly selected boundary model.

## Where the custom implementation may be

A read-only public GitHub inspection on 2026-09-10 found only `main` and `RyanZhu` branches in `marsbarslars/WarpX-Anisotropy`, with no C++ source or patch implementing this boundary. Lars's public repository list did not include a WarpX source fork. The `chore/uv-project-scaffold` branch of `marsbarslars/Plasma-Hackathon` points its `vendor/warpx` submodule at upstream WarpX `development`.

The shared file `simulations/Mirror/amrex_source_code_navigation.sh` points to `~/miniconda3/pkgs/amrex-26.06-nompi_dp_hec1556b_101/include`. This is evidence that Lars was navigating an installed AMReX package, but it does not prove where his edits live. His changes may be local or otherwise not publicly accessible. No branch or package has been modified as part of this audit.

The useful item to obtain from Lars is the actual source patch or repository/branch/commit implementing `Periodic_ReflectParVel`, together with the matching WarpX and AMReX versions and a description of the intended mapping. The name alone does not establish whether it flips a Cartesian component, reflects velocity relative to the local magnetic field, rotates coordinates between ends, or changes particle speed.

## Boundary model for this exploratory workflow

The new workflow explicitly selects ordinary periodic boundaries in z and preserves the original custom input separately. Standard periodic wrapping preserves the entire velocity vector at the crossing. It does not implement Lars's proposed direction or speed change.

The available magnetic-field audit found a sign change in transverse field components between the two z ends. That mismatch matters for repeated passage across a periodic seam; it is not evidence of a correctly joined magnetic geometry. The planned kinetic startup uses an interior background region, `0.25 < z < 4.75 m` and `r < 0.75 m`, and a short 10 ns horizon to examine initialization, electron/ion dynamics, and collisions. Starting away from the seam reduces exposure during that short test; it does not repair the seam or validate long-time confinement.

A velocity direction change at fixed speed preserves kinetic energy. A speed change changes kinetic energy and needs a specified physical or coordinate mapping; it must not be introduced as an undocumented periodic-boundary tweak. Boundary-crossing transformations and long-time loss/confinement conclusions remain pending Lars's implementation and its validation. The absorbing transverse/embedded boundaries and particle injection also mean that the domain is not generally a closed particle or energy system.

## References

- [Shared input at the inspected commit](https://github.com/marsbarslars/WarpX-Anisotropy/blob/b8acce9b91d65b0ddc5102550e8e29c43ffacc02/simulations/Mirror/inputs_mirror.txt)
- [WarpX boundary enum at the installed commit](https://github.com/BLAST-WarpX/warpx/blob/bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc/Source/Utils/WarpXAlgorithmSelection.H#L144)
- [WarpX particle-boundary parser at the installed commit](https://github.com/BLAST-WarpX/warpx/blob/bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc/Source/Particles/ParticleBoundaries.cpp#L53)
- [AMReX periodic position mapping at the installed commit](https://github.com/AMReX-Codes/amrex/blob/057940244648b82908cfc486f07ec796bba2b07f/Src/Particle/AMReX_ParticleUtil.H#L423)
