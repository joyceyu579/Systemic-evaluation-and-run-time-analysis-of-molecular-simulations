import yu_joyce_monte_carlo as mc

# Set simulation parameters
reduced_temperature = 0.9
num_steps = 5000
max_displacement = 0.1 
cutoff = 3
freq = 1000

# Read initial coordinates
coordinates, box_length = mc.read_xyz("../data/sample_config1.txt")
#coordinates, box_length = mc.generate_cubic_lattice_config(800, 1)

mc.run_simulation(coordinates, box_length, reduced_temperature, cutoff, num_steps, max_displacement)