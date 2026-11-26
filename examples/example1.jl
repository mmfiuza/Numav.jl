using Numav

mesh_path = "example1.bdf"
result_path = "result.h5"

# Create the simulation object
s = create_simulation(
    numerical_method = Fem,
    equation = Helmholtz,
    element_shape = Tetrahedron,
    element_order = Linear
)

# Set simulated frequencies
set_frequency!(s, min=20, max=200, step=1)

# Load mesh
load_mesh!(s, mesh_path)

# Define air
rho0 = 1.2
c0= 343
add_volume_material!(s, region=1, density=rho0, speed_of_sound=c0)

# Define material as an equivalent fluid
rho1 = 1 - 1im
c1 = 400 + 400im
add_volume_material!(s, region=2, density=rho1, speed_of_sound=c1)

# Add a surface impedance boundary condition
Zs = 4000 + 4000im
add_surface_material!(s, region=11, specific_acoustic_impedance=Zs)

# Add a monopole volume velocity source at a point
Q(f) = 0.01/f
add_sound_source!(s, coordinates=[0.0, 5.0, 0.0], volume_velocity=Q)

# Add a particle velocity source on a surface
U(f) = 0.015/f
add_sound_source!(s, region=10, particle_velocity=U)

# Add pressure source
P(f) = 0.002f
add_sound_source!(s, coordinates=[4.0, 5.0, 0.0], pressure=P)

# Set output file and run
set_result_export_path!(s, result_path)
run!(s)

# Plot result
plot_pressure_field(result_path)
