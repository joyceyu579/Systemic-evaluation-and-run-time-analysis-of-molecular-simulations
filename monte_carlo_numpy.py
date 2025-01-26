"""
NumPy Implementation of MC simulation code.
"""

import math
import random
import numpy as np

def calculate_LJ_np(r_ij):
    r6_term = np.power(1./r_ij, 6.)
    r12_term = np.power(r6_term, 2.)
    pairwise_energy = 4 * (r12_term - r6_term)
    return pairwise_energy

def calculate_distance_np(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.
    Parameters
    ----------
    coord1, coord2: np.ndarray
        The atomic coordinates

    Returns
    -------
    distance: float
        The distance between the two points.
    """

    dim_dist = coord2 - coord1
    
    if box_length:
        dim_dist = dim_dist - box_length * np.round(dim_dist / box_length)

    if dim_dist.ndim < 2:
        dim_dist = dim_dist.reshape(1, -1)

    dim_dist_sq = dim_dist ** 2

    dim_dist_sum = dim_dist_sq.sum(axis=1)

    distance = np.sqrt(dim_dist_sum)

    return distance
    
def calculate_total_energy(coordinates, box_length, cutoff):
    """
    Calculate the total Lennard Jones energy of a system of particles.

    Parameters
    ----------
    coordinates : list
        Nested list containing particle coordinates.

    Returns
    -------
    total_energy : float
        The total pairwise Lennard Jones energy of the system of particles.
    """

    total_energy = 0

    num_atoms = len(coordinates)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):

            dist_ij = calculate_distance_np(
                coordinates[i], coordinates[j], box_length=box_length
            )

            if dist_ij < cutoff:
                interaction_energy = calculate_LJ_np(dist_ij)
                total_energy += interaction_energy

    return total_energy

def read_xyz(filepath):
    """
    Reads coordinates from an xyz file.
    
    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.
       
    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates

    """
    
    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())
        coordinates = f.readlines()
    
    atomic_coordinates = []
    
    for atom in coordinates:
        split_atoms = atom.split()
        
        float_coords = []
        
        # We split this way to get rid of the atom label.
        for coord in split_atoms[1:]:
            float_coords.append(float(coord))
            
        atomic_coordinates.append(float_coords)
        
    
    return atomic_coordinates, box_length

def calculate_tail_correction(num_particles, box_length, cutoff):
    """
    Calculate the long range tail correction
    """

    const1 = (8 * math.pi * num_particles**2) / (3 * box_length**3)
    const2 = (1 / 3) * (1 / cutoff) ** 9 - (1 / cutoff) ** 3

    return const1 * const2

def generate_cubic_lattice(num_atoms, density):

    # Calculate box length based on number of atoms and density. 
    volume = num_atoms / density
    box_length = math.pow(volume, (1./3.))

    #Calculate the upper bound of cube size. 
    # Our approach will be to place atoms until 
    # we place all needed. For this, we need 
    # to determine the maximum number of atoms on each side.
    max_side = math.ceil(math.pow(num_atoms, (1./3.)))

    # Determine spacing based on number of atoms and box length.
    spacing = box_length / max_side # units length / atom

    coordinates = []
    count = 0

    for i in range(max_side):
        for j in range(max_side):
            for k in range(max_side):
                coordinates.append([i*spacing, j*spacing, k*spacing])
                count += 1
                if count == num_atoms:
                    return coordinates, box_length

def accept_or_reject(delta_e, beta):
    """
    Accept or reject based on an energy and beta (inverse temperature) measurement.

    Parameters
    ----------
    delta_e : float
        An energy change in reduced units.
    beta : float
        Inverse temperature.

    Returns
    -------
    bool
        True to accept move, False to reject
    """

    if delta_e <= 0.0:
        accept = True
    else: 
        random_number = random.random()
        p_acc = math.exp(-delta_e*beta)

        if random_number < p_acc:
            accept = True
        else:
            accept = False

    return accept

def calculate_pair_energy_np(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system).

    Parameters
    ----------
    coordinates : np.ndarray
        The coordinates for all particles in the system.
    i_particle : int
        The particle index for which to calculate the energy.
    cutoff : float
        The simulation cutoff. Beyond this distance, the interactions are not calculated. 

    Returns
    -------
    e_total : float
        The pairwise interaction energy of the i_th particle with all other particles.
    """
    
    i_position = coordinates[i_particle]

    j_positions = coordinates

    r_ij = calculate_distance_np(i_position, j_positions, box_length)
    
    less_than_cutoff = r_ij[np.where(r_ij < cutoff) and np.where(r_ij != 0)]

    e_pair = calculate_LJ_np(less_than_cutoff)
    e_total = e_pair.sum(axis=0)

    return e_total

# Create a function called run_simulation that takes in simulation parameters and runs the simulation.
def run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000):

    # Calculated quantities
    beta = 1 / reduced_temperature
    num_particles = len(coordinates)

    delta_e = 0

    total_energy = calculate_total_energy(coordinates, box_length, cutoff)
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff)

    steps = []
    energies = []

    for step in range(num_steps):

        # 1. Randomly pick one of N particles
        random_particle = random.randrange(num_particles)

        # 2. Calculate the interaction energy of the selected particle with the system and store this value.
        current_energy = calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)

        # 3. Generate a random x, y, z displacement
        x_rand = random.uniform(-max_displacement, max_displacement)
        y_rand = random.uniform(-max_displacement, max_displacement)
        z_rand = random.uniform(-max_displacement, max_displacement)

        # 4. Modify the coordinate of Nth particle by the generated displacements.
        coordinates[random_particle][0] += x_rand
        coordinates[random_particle][1] += y_rand
        coordinates[random_particle][2] += z_rand

        # 5. Calculate the interaction energy of the selected particle with the system and store this value (using the new position).
        proposed_energy = calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)

        # 6. Calculate the change in energy and decide to accept or reject the change. 
        delta_e = proposed_energy - current_energy

        accept = accept_or_reject(delta_e, beta)

        # 7. If accept, keep particle movement. If reject, move particle back.
        if accept:
            total_energy += delta_e
        else: 
            # Move is not accepted, roll back coordinates
            coordinates[random_particle][0] -= x_rand
            coordinates[random_particle][1] -= y_rand
            coordinates[random_particle][2] -= z_rand

        # 8. Print the energy if the step is a multiple of freq
        if step % freq == 0:
            print(step, total_energy/num_particles)
            steps.append(step)
            energies.append(total_energy)


