import openpmd_api as io

series = io.Series(
    "diags/diag2/openpmd_000005.bp5",
    io.Access.read_only
)

for it in series.iterations:
    print("Iteration:", it)

    iteration = series.iterations[it]

    print("Meshes:")
    for name in iteration.meshes:
        print("  ", name)

    print("Particles:")
    for name in iteration.particles:
        print("  ", name)

    if "protons" in iteration.particles:
        protons = iteration.particles["protons"]

        print("\nProton records:")
        for name in protons:
            print("  ", name)

series = io.Series(
    "diags/diag2/openpmd_000005.bp5",
    io.Access.read_only
)

iteration = series.iterations[5]
protons = iteration.particles["protons"]

for record_name in protons:
    record = protons[record_name]
    print(f"\n{record_name}:")
    for component_name in record:
        component = record[component_name]
        print(
            f"  {component_name}: "
            f"shape={component.shape}, "
            f"dtype={component.dtype}"
        )

series.close()