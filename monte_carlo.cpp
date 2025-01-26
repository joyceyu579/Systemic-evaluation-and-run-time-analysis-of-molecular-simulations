#include <iostream>
#include <math.h>
#include <fstream>
#include <array>
#include <vector>
#include <utility>
#include <random>
#include <chrono>


// A Global! Probably shouldn't be used in real code
std::default_random_engine re;

/*! Generate a random double within a given range */
double random_double(double lower_bound, double upper_bound)
{
   std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
   return dist(re);
}

/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
int random_integer(int lower_bound, int upper_bound)
{           
   //dist will return [a,b] but we want [a,b)
   std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
   return dist(re);
}  

// Make some types more convenient
typedef std::array<double, 3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;

std::pair<Coordinates, double> read_xyz(std::string file_path)
{
    // Opens up a file stream for input
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if(!infile.is_open())
    {   
        throw std::runtime_error("File path in read_xyz does not exist!");
    }
    
    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;
    
    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;
    
    // now the number of atoms
    infile >> num_atoms;
    
    // Uncomment to help troubleshoot
    //std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;
    
    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;
    
    for(int i = 0; i < num_atoms; i++)
    {   
        AtomCoord coord;
        
        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];
        
        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}


double calculate_LJ(double r_ij)
{
    double r6_term = pow((1/r_ij), 6.);
    double r12_term = pow(r6_term, 2.);
    double pairwise_energy = 4 * (r12_term - r6_term);

    return pairwise_energy;
}

double calculate_tail_correction(int num_particles, double box_length, double cutoff)
{
    double const1 = (8 * M_PI * pow(num_particles, 2.)) / (3 * pow(box_length, 3.));
    double const2 = (1. / 3.) * pow((1. / cutoff), 9.) - pow((1. / cutoff), 3.);

    return const1 * const2;
}

bool accept_or_reject(double delta_e, double beta)
{
    bool accept;

    if(delta_e <= 0)
    {
        accept = true;
    }
    else
    {
        double random_number = random_double(0.0, 1.0);
        double p_acc = exp(-delta_e*beta);

        if(random_number < p_acc)
        {
            accept = true;
        }
        else
        {
            accept = false;
        }
    }

    return accept;
}

double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length)
{
    double distance = 0;
    double dim_dist;

    for(int i = 0; i < 3; i++)
    {
        dim_dist = coord1[i] - coord2[i];

        if(box_length)
        {
            dim_dist = dim_dist - box_length * round(dim_dist / box_length);
        }

        dim_dist = pow(dim_dist, 2.);
        distance += dim_dist;
    }

    distance = sqrt(distance);
    return distance;
}

std::pair<Coordinates, double> generate_rand_config(int num_particles, double density)
{
    double box_length = pow(num_particles/density, 1/3.);
    Coordinates atomic_coordinates;

    for(int i = 0; i < num_particles; i++)
    {
        double x = random_double(-box_length, box_length);
        double y = random_double(-box_length, box_length);
        double z = random_double(-box_length, box_length);

        AtomCoord coords = {x, y, z};
        atomic_coordinates.push_back(coords); 
    }

    return std::make_pair(atomic_coordinates, box_length);
}

std::pair<Coordinates, double> generate_cubic_lattice(int num_atoms, double density)
{
    double volume = num_atoms / density;
    double box_length = pow(volume, 1./3.);

    int max_side = ceil(pow(num_atoms, (1./3.)));

    double spacing = box_length / max_side;

    Coordinates coordinates;
    int count = 0;

    for(int i = 0; i < max_side; i++)
    {
        for(int j = 0; j < max_side; j++)
        {
            for(int k = 0; k < max_side; k++)
            {
                coordinates.push_back({i * spacing, j * spacing, k * spacing});
                count++; 
                if(count == num_atoms)
                {
                    break;
                }
            }
        }
    }

    return std::make_pair(coordinates, box_length);
}

double calculate_pair_energy(Coordinates & coordinates, double i_particle, double box_length, double cutoff)
{
    double e_total = 0.0;
    AtomCoord i_position = coordinates[i_particle];

    double num_atoms = coordinates.size();

    for(int j_particle = 0; j_particle < num_atoms; j_particle++)
    {
        if(i_particle != j_particle)
        {
            AtomCoord j_position = coordinates[j_particle];

            double r_ij = calculate_distance(i_position, j_position, box_length);

            if(r_ij < cutoff)
            {
                double e_pair = calculate_LJ(r_ij);
                e_total += e_pair;
            }
        }
    }

    return e_total;
}

double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff)
{
    double total_energy = 0.0;

    double num_atoms = coordinates.size();

    for(int i = 0; i < num_atoms; i++)
    {
        for(int j = i+1; j < num_atoms; j++)
        {
            double dist_ij = calculate_distance(
                coordinates[i], coordinates[j], box_length=box_length
            );

            if(dist_ij < cutoff)
            {
                double interaction_energy = calculate_LJ(dist_ij);
                total_energy += interaction_energy;
            } 
        }
    }

    return total_energy;
}

std::pair<std::vector<int>, std::vector<double>> run_simulation(Coordinates coordinates, double box_length, double cutoff, double reduced_temperature, int num_steps, double max_displacement=0.1, int freq=1000)
{
    double beta = 1 / reduced_temperature;
    int num_particles = coordinates.size();

    double delta_e = 0.0;

    double total_energy = calculate_total_energy(coordinates, box_length, cutoff);
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);

    std::vector<int> steps = {};
    std::vector<double> energies = {};

    for(int step = 0; step < num_steps; step++)
    {
        int random_particle = random_integer(0, num_particles);

        double current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;

        double proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        delta_e = proposed_energy - current_energy;

        bool accept = accept_or_reject(delta_e, beta);

        if(accept)
        {
            total_energy += delta_e;
        }
        else
        {
            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
        }

        if(step % freq == 0)
        {
            std::cout << step << " " << total_energy/num_particles << std::endl;
            steps.push_back(step);
            energies.push_back(total_energy);
        }
    }
    return std::make_pair(steps, energies);
}

int main(void)
{
    
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());
    
    std::pair<Coordinates, double> xyz_info = read_xyz("../data/sample_config1.txt");

    Coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;
    
    
    double reduced_temperature = 0.9;
    int num_steps = 100000;
    double cutoff = 3.0;
    
    run_simulation(coords, box_length, cutoff, reduced_temperature, num_steps);

   return 0;
}