// Include these at the top of your file
#include <iostream> // for std::cout, std::endl
using namespace std;
#include <fstream> // for reading/writing files
#include <array>   // for std::array
#include <vector>  // for std::vector
#include <utility> // for std::pair

typedef array<double, 3> AtomCoord; 
typedef vector<AtomCoord> coordinates;

#include <random> // for random numbers
#include <chrono> // for generating the seed


std::pair<coordinates, double> read_xyz(std::string file_path)
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
    coordinates coords;
    
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


double calculate_LJ(float r_ij)
{
    float r6_term;
    float r12_term;
    float pairwise_energy;

    r6_term = pow(1/r_ij, 6);
    r12_term = pow(r6_term, 2);
    pairwise_energy = 4* (r12_term - r6_term);

    return pairwise_energy;
}

double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length)
{
    double distance = 0; 
    double dim_dist = 0;

    for(int i = 0; i < 3; i++)
    {
        dim_dist = coord1[i] - coord2[i];
        if (box_length)
        {
            dim_dist = dim_dist - box_length * round(dim_dist / box_length); 
        }
        dim_dist = pow(dim_dist,2);
        distance += dim_dist;
    }
    distance = sqrt(distance);
    return distance;
}

double calculate_total_energy(coordinates MyCoords, double box_length, double cutoff)

{
    double total_energy = 0;
    double num_atoms = MyCoords.size();

    for(double i = 0; i < num_atoms; i++)
    {
        for(double j = i+1; j < num_atoms; j++)
        {
            double dist_ij = calculate_distance(MyCoords[i], MyCoords[j], box_length = box_length);
            if(dist_ij < cutoff)
            {
                double interaction_energy = calculate_LJ(dist_ij);
                total_energy += interaction_energy;
            }
    }}

    return total_energy; 
}

double calculate_tail_correction(int num_particles, float box_length, float cutoff)
{
    float const1;
    float const2;

    const1 = (8 * 3.1415 /*pi*/ * pow(num_particles,2) / (3 * pow(box_length, 3)));
    const2 = (1/3) * pow((1 / cutoff), 9) - pow((1/ cutoff), 3);

    return const1 * const2;
}

pair <coordinates, double> generate_cubic_lattice(double num_particles, double density)
{
    double box_length = pow(num_particles/density, 1/3.); // need period at end. 
    coordinates atomic_coordinates; 
    double max_side = ceil(pow(num_particles, (1./3.)));
    double spacing = box_length / max_side;
    int count = 0;

    for(int i = 0; i < max_side; i++)
    {
        for (int j = 0; j < max_side; j++)
        {
            for (int k = 0; k < max_side; k++)
            {
                AtomCoord coords = {i*spacing, j*spacing, k*spacing};
                atomic_coordinates.push_back(coords);
                count += 1;
                if(count == num_particles)
                {
                    break;
                }
            }
        }
    }
    return make_pair(atomic_coordinates, box_length); 
}

pair <coordinates, double> generate_rand_config(double num_particles, double density)
{
    double box_length = pow(num_particles/density, 1/3.); // need period at end. 
    coordinates atomic_coordinates; 

    for(int i = 0; i < num_particles; i++)
    {
        double x = random_double(-box_length, box_length);
        double y = random_double(-box_length, box_length);
        double z = random_double(-box_length, box_length);

        AtomCoord coords = {x, y, z};
        atomic_coordinates.push_back(coords);
    }
    return make_pair(atomic_coordinates, box_length); 
}

bool accept_or_reject(float delta_e, float beta)
{
    bool accept;
    double random_number;
    float p_acc;

    if(delta_e <= 0.0)
    {
        accept = true;
    } else { 
        random_number = random_double(0, 1);
        p_acc = exp(-delta_e * beta);

        if(random_number < p_acc)
        {
            accept = true;
        }else{
            accept = false;
        }
    }
    return accept;
}

double pair_energy(coordinates MyCoords, double i_particle, double box_length, double cutoff)
{
    double e_total = 0.0;
    AtomCoord i_position = MyCoords[i_particle];
    double num_atoms = MyCoords.size();
    for(double j_particle; j_particle < num_atoms; j_particle++)
    {
        if(i_particle != j_particle)
        {
            AtomCoord j_position = MyCoords[j_particle];
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

pair < vector<int>, vector<double>>run_simulation(coordinates MyCoords, double box_length, double cutoff, double reduced_temperature, int num_steps, double max_displacement=0.1, int freq = 1000)
{
    double beta = 1 / reduced_temperature;
    double num_particles = MyCoords.size();
    
    double total_energy = calculate_total_energy(MyCoords, box_length, cutoff);
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);

    vector<int> steps = {};
    vector<double> energies = {}; 
    double delta_e = 0;
    
    for(int step = 0; step < num_steps; step++)
    {
        double random_particle = random_double(0, num_particles);
        double current_energy = pair_energy(MyCoords, random_particle, box_length, cutoff);
        //#3 Generate random particle
        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);
        //#4 Modify the coordinate of nth particle by the generated displacements. 
        MyCoords[random_particle][0] += x_rand;
        MyCoords[random_particle][1] += y_rand;
        MyCoords[random_particle][2] += z_rand;
        //#5 
        float proposed_energy = pair_energy(MyCoords, random_particle, box_length, cutoff);
        //#6
        float delta_e = proposed_energy - current_energy;
        float accept = accept_or_reject(delta_e, beta);
        //#7
        if(accept)
        {
            total_energy += delta_e;
        }else{
            MyCoords[random_particle][0] -= x_rand;
            MyCoords[random_particle][1] -= y_rand;
            MyCoords[random_particle][2] -= z_rand;
        }
        if (step % freq == 0)
        {
            cout << "steps: " << step  << "energy: " << total_energy/num_particles << endl;
            steps.push_back(step);
            energies.push_back(total_energy);
        }
    }
return make_pair(steps, energies);
}


int main()
{
    // Initialize random number generation based on time
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());

    std::pair<coordinates, double> xyz_info = read_xyz("../data/sample_config1.txt");

    coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;

    run_simulation(coords, box_length, 3.0, 0.9, 5000);
}