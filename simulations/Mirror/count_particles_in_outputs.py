import glob
import openpmd_api as io

files = sorted(glob.glob("diags/diag1/openpmd_*.bp5"))

for filename in files:
    series = io.Series(filename, io.Access.read_only)

    for iteration, it in series.iterations.items():

        print(f"{filename}: iteration {iteration}")

        for species_name, species in it.particles.items():

            # Use position/x to determine the number of particles
            x = species["position"]["x"]

            n_particles = x.shape[0]

            print(f"    {species_name}: {n_particles} particles")

    series.close()