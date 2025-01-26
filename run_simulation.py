import monte_carlo as mc

# Set simulation parameters
reduced_temperature = 0.9 
num_steps = 50000
cutoff = 3

# Read or generate initial coordinates
#coordinates, box_length = mc.read_xyz("../data/sample_config1.txt")
#coordinates, box_length = mc.generate_cubic_lattice(800, 0.8)

config_file1 = "../data/sample_config1.txt"

coordinates, box_length = mc.read_xyz(config_file1)

# Run simulation
mc.run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps)