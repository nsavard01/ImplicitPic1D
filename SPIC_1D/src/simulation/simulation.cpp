#include "simulation/simulation.hpp"
#include "globals/plasma_functions.hpp"


simulation::simulation() {
    // Get scheme type from the file
    int scheme_type = 0;
    if (mpi_vars::mpi_rank == 0) {
        std::ifstream file("../inputs/initial_setup.inp");
        if (!file) {
            std::cout << "Error: Unable to open file initial_setup.inp" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string line;
        std::getline(file, line);
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> scheme_type;
        file.close();
        if (scheme_type == 0) {
            std::cout << "Scheme type: 0 (MC-PIC)" << std::endl;
        } else if (scheme_type == 1) {
            std::cout << "Scheme type: 0 (EC-PIC)" << std::endl;
        } else if (scheme_type == 2) {
            std::cout << "Scheme type: 2 (I-NGP)" << std::endl;
        } else if (scheme_type == 3) {
            std::cout << "Scheme type: 3 (I-CIC)" << std::endl;
        } else {
            std::cerr << "Error: Unknown scheme type." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    MPI_Bcast(&scheme_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
    this->scheme_type = scheme_type;
    initialize_pcg(false); // Initialize the PCG RNG with a non-deterministic seed
    // Generate objects serially except when needed
    this->world = create_domain_from_file("../inputs/geometry.inp", scheme_type);
    this->world->print_out();
    this->charged_particle_list = read_charged_particle_inputs("../inputs/charged_particles/", *this->world); 
    this->del_t = 0.2 / get_plasma_frequency(this->charged_particle_list[0].average_temperature, this->charged_particle_list[0].average_density); // Time step
    this->field_solver = read_voltage_inputs("../inputs/geometry.inp", this->scheme_type, *this->world);
    // create openmp parallel, with thread id passed around
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        this->field_solver->deposit_charge_density(this->charged_particle_list, thread_id);
        #pragma omp master
        {
            this->field_solver->solve_potential(0.0, *this->world);
            this->field_solver->make_EField(*this->world);
        }
        #pragma omp barrier
        this->field_solver->push_particles(thread_id, this->del_t, this->charged_particle_list, *this->world);
    }
}
