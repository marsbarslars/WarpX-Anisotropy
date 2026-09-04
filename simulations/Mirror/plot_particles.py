import numpy as np
import pyvista as pv
from openpmd_viewer import OpenPMDTimeSeries


# ============================================================
# Settings
# ============================================================

diag_dir = "diags/diag1"
output_file = "movies/periodic.mp4"

fps = 30

particle_size = 6
particle_color = "red"

# Number of iterations kept in each particle trail
trail_length = 10

trail_width = 2.0
trail_color = "red"

# Number of magnetic field lines
n_field_lines = 20

# Length of field lines
streamline_length = 10.0

streamline_max_steps = 5000
streamline_terminal_speed = 1e-8

bfield_line_width = 1.5
bfield_color = "white"


# ============================================================
# Load OpenPMD data
# ============================================================

series = OpenPMDTimeSeries(diag_dir)

iterations = series.iterations

print(
    f"Found {len(iterations)} iterations"
)


# ============================================================
# Load magnetic field and particles
# ============================================================

B_data = []
particle_data = []

x_coords = None
y_coords = None
z_coords = None


for i, it in enumerate(iterations):

    print(
        f"Loading iteration {it} "
        f"({i + 1}/{len(iterations)})"
    )

    # --------------------------------------------------------
    # Magnetic field
    # --------------------------------------------------------

    Bx, info = series.get_field(
        "B",
        coord="x",
        iteration=it
    )

    By, _ = series.get_field(
        "B",
        coord="y",
        iteration=it
    )

    Bz, _ = series.get_field(
        "B",
        coord="z",
        iteration=it
    )


    if x_coords is None:

        x_coords = np.asarray(info.x)
        y_coords = np.asarray(info.y)
        z_coords = np.asarray(info.z)


    B_data.append(
        np.stack(
            (
                Bx,
                By,
                Bz
            ),
            axis=0
        )
    )


    # --------------------------------------------------------
    # Particle data
    # --------------------------------------------------------

    try:

        x, y, z, ids = series.get_particle(
            ["x", "y", "z", "id"],
            iteration=it,
            species="protons"
        )

        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        ids = np.asarray(ids)

    except Exception:

        x = np.empty(0)
        y = np.empty(0)
        z = np.empty(0)

        ids = np.empty(
            0,
            dtype=np.int64
        )


    particle_data.append(
        {
            "id": ids,
            "x": x,
            "y": y,
            "z": z
        }
    )


B_data = np.asarray(B_data)


print("Finished loading data")


# ============================================================
# Find unique particle IDs
# ============================================================

all_particle_ids = set()

for data in particle_data:

    all_particle_ids.update(
        data["id"].tolist()
    )


print(
    f"Found {len(all_particle_ids)} "
    f"unique particle IDs"
)


# ============================================================
# Create B-field grid
# ============================================================

X, Y, Z = np.meshgrid(
    x_coords,
    y_coords,
    z_coords,
    indexing="ij"
)


nx = len(x_coords)
ny = len(y_coords)
nz = len(z_coords)


grid = pv.StructuredGrid()


grid.points = np.column_stack(
    (
        X.ravel(),
        Y.ravel(),
        Z.ravel()
    )
)


grid.dimensions = (
    nx,
    ny,
    nz
)


grid["B"] = np.column_stack(
    (
        B_data[0, 0].ravel(),
        B_data[0, 1].ravel(),
        B_data[0, 2].ravel()
    )
)


# ============================================================
# Create candidate field-line seeds
#
# All seeds lie in the XZ plane:
#
#     y = 0
#
# We generate many candidates and keep the first 10
# that actually produce a streamline.
# ============================================================

n_candidate_x = 20
n_candidate_z = 20


candidate_x = np.linspace(
    x_coords.min(),
    x_coords.max(),
    n_candidate_x
)


candidate_z = np.linspace(
    z_coords.min(),
    z_coords.max(),
    n_candidate_z
)


candidate_seeds = []


for x in candidate_x:

    for z in candidate_z:

        candidate_seeds.append(
            [
                x,
                0.0,
                z
            ]
        )


candidate_seeds = np.asarray(
    candidate_seeds
)


print(
    "Generating magnetic field lines..."
)


# ============================================================
# Find 10 usable field lines
# ============================================================

valid_streamlines = []

used_seeds = []


for seed in candidate_seeds:

    if len(valid_streamlines) >= n_field_lines:
        break


    seed_poly = pv.PolyData(
        np.asarray(
            [seed]
        )
    )


    try:

        lines = (
            grid.streamlines_from_source(
                seed_poly,
                vectors="B",
                integration_direction="both",
                max_length=streamline_length,
                max_steps=streamline_max_steps,
                terminal_speed=streamline_terminal_speed
            )
        )

    except Exception:

        continue


    if lines.n_points > 1:

        # ----------------------------------------------------
        # Project streamline onto XZ plane
        # ----------------------------------------------------

        points = lines.points.copy()

        points[:, 1] = 0.0

        lines.points = points


        valid_streamlines.append(
            lines
        )

        used_seeds.append(
            seed
        )


print(
    f"Found {len(valid_streamlines)} "
    f"usable magnetic field lines"
)


# ============================================================
# Combine field lines
# ============================================================

if len(valid_streamlines) > 0:

    field_line_mesh = (
        pv.merge(
            valid_streamlines
        )
    )

else:

    field_line_mesh = None


# ============================================================
# Create plotter
# ============================================================

plotter = pv.Plotter(
    window_size=(1200, 800)
)


plotter.set_background(
    "black"
)


# ============================================================
# Add magnetic field lines
# ============================================================

if field_line_mesh is not None:

    bfield_actor = plotter.add_mesh(
        field_line_mesh,
        color=bfield_color,
        line_width=bfield_line_width,
        name="bfield"
    )


# ============================================================
# Particle mesh
# ============================================================

particle_mesh = pv.PolyData(
    np.array(
        [
            [0.0, 0.0, 0.0]
        ]
    )
)


particle_actor = plotter.add_mesh(
    particle_mesh,
    color=particle_color,
    point_size=particle_size,
    render_points_as_spheres=True,
    name="particles"
)


# ============================================================
# Trail mesh
# ============================================================

trail_mesh = pv.PolyData(
    np.array(
        [
            [0.0, 0.0, 0.0]
        ]
    )
)


trail_actor = plotter.add_mesh(
    trail_mesh,
    color=trail_color,
    line_width=trail_width,
    name="trails"
)


# ============================================================
# XZ camera
# ============================================================

xmin = x_coords.min()
xmax = x_coords.max()

zmin = z_coords.min()
zmax = z_coords.max()


xcenter = 0.5 * (
    xmin + xmax
)

zcenter = 0.5 * (
    zmin + zmax
)


domain_size = max(
    xmax - xmin,
    zmax - zmin
)


# Look directly along Y
plotter.camera.position = (
    xcenter,
    10.0 * domain_size,
    zcenter
)


plotter.camera.focal_point = (
    xcenter,
    0.0,
    zcenter
)


plotter.camera.up = (
    0,
    0,
    1
)


plotter.camera.parallel_projection = True


plotter.show_axes()


# ============================================================
# Text
# ============================================================

plotter.add_text(
    f"Iteration {iterations[0]}",
    position="upper_left",
    font_size=16,
    name="iteration_text"
)


# ============================================================
# Particle history
# ============================================================

history = {}


# ============================================================
# Start MP4
# ============================================================

print(
    "Writing MP4..."
)


plotter.open_movie(
    output_file,
    framerate=fps
)


# ============================================================
# Animation loop
# ============================================================

for frame, it in enumerate(iterations):

    print(
        f"Rendering frame "
        f"{frame + 1}/{len(iterations)} "
        f"(iteration {it})"
    )


    data = particle_data[
        frame
    ]


    ids = data["id"]
    xs = data["x"]
    zs = data["z"]


    # ========================================================
    # Update particle history
    # ========================================================

    for pid, x, z in zip(
        ids,
        xs,
        zs
    ):

        position = np.array(
            [
                x,
                z
            ]
        )


        if pid not in history:

            history[pid] = []


        history[pid].append(
            position
        )


        # Keep last 10 observations
        if len(history[pid]) > trail_length:

            history[pid] = (
                history[pid][
                    -trail_length:
                ]
            )


    # ========================================================
    # Current particle positions
    # ========================================================

    if len(ids) > 0:

        current_points = np.column_stack(
            (
                xs,
                np.zeros_like(xs),
                zs
            )
        )


        particle_mesh.points = (
            current_points
        )

    else:

        particle_mesh.points = np.array(
            [
                [
                    1e30,
                    1e30,
                    1e30
                ]
            ]
        )


    # ========================================================
    # Build combined trail mesh
    # ========================================================

    trail_points = []
    trail_lines = []

    point_counter = 0


    for positions in history.values():

        if len(positions) < 2:
            continue


        positions = np.asarray(
            positions
        )


        for j in range(
            len(positions) - 1
        ):

            x0, z0 = positions[j]

            x1, z1 = positions[j + 1]


            trail_points.append(
                [
                    x0,
                    0.0,
                    z0
                ]
            )


            trail_points.append(
                [
                    x1,
                    0.0,
                    z1
                ]
            )


            trail_lines.extend(
                [
                    2,
                    point_counter,
                    point_counter + 1
                ]
            )


            point_counter += 2


    # ========================================================
    # Update trail mesh
    # ========================================================

    if len(trail_points) > 0:

        new_trail_mesh = pv.PolyData(
            np.asarray(
                trail_points
            ),
            lines=np.asarray(
                trail_lines,
                dtype=np.int64
            )
        )


        trail_actor.mapper.dataset = (
            new_trail_mesh
        )


    # ========================================================
    # Update text
    # ========================================================

    plotter.add_text(
        f"Iteration {it}",
        position="upper_left",
        font_size=16,
        name="iteration_text"
    )


    # ========================================================
    # Write frame
    # ========================================================

    plotter.write_frame()


# ============================================================
# Finish
# ============================================================

plotter.close()


print()
print(
    f"Saved: {output_file}"
)