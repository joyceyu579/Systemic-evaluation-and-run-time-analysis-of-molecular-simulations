import monte_carlo_numpy as mc
import numpy as np

# Set simulation parameters
reduced_temperature = 0.9 
num_steps = 100000
cutoff = 3

# Read or generate initial coordinates
coordinates, box_length = mc.read_xyz("../data/sample_config1.txt")

coordinates_np = np.array(coordinates)

#print(len(coordinates_np))
mc.run_simulation(coordinates_np, box_length, cutoff, reduced_temperature, num_steps)