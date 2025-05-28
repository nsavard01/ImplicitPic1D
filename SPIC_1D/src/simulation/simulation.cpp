#include "simulation/simulation.hpp"
#include "globals/plasma_functions.hpp"


simulation::simulation() {
    // Get scheme type from the file
    int max_threads = omp_get_max_threads();

    // Get number of OpenMP threads from the file
    int number_omp_threads;
    if (mpi_vars::mpi_rank == 0) {
        std::ifstream file("../inputs/initial_setup.inp");
        if (!file) {
            std::cout << "Error: Unable to open file initial_setup.inp" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string line;
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> number_omp_threads;
        file.close();
        if (number_omp_threads > max_threads) {
            std::cerr << "Error: Requested number of OpenMP threads (" << number_omp_threads
                    << ") exceeds the maximum available (" << max_threads << ")." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&number_omp_threads, 1, MPI_INT, 0, MPI_COMM_WORLD);
    omp_set_num_threads(number_omp_threads); // Set the number of OpenMP threads
    for (int i = 0; i < mpi_vars::mpi_size; i++) {
        if (i == mpi_vars::mpi_rank) {
            std::cout << "MPI rank: " << mpi_vars::mpi_rank << " using " << omp_get_max_threads() << " OpenMP threads out of maximum of " << max_threads << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
}


void simulation::setup() {
    // create openmp parallel, with thread id passed around
    int scheme_type = 0;
    double del_t_temp, del_t_fraction;
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
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> del_t_fraction >> del_t_temp;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> this->simulation_time;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> this->averaging_time;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> this->number_diagnostics;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> this->save_file_folder;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> this->save_file_path;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        std::string rest;
        iss >> rest;
        this->restarted_simulation = (rest == "yes" || rest == "Yes");
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
    MPI_Bcast(&del_t_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&del_t_fraction, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&this->simulation_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&this->averaging_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&this->number_diagnostics, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&this->restarted_simulation, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    this->scheme_type = scheme_type;
    initialize_pcg(false); // Initialize the PCG RNG with a non-deterministic seed
    // Generate objects serially except when needed
    this->world = create_domain_from_file("../inputs/geometry.inp", scheme_type);
    this->world->print_out();
    this->charged_particle_list = read_charged_particle_inputs("../inputs/charged_particles/", *this->world); 
    this->target_particle_list = read_target_particle_inputs("../inputs/target_particles/", *this->world);
    this->null_collider_list = read_null_collision_inputs("../inputs/collisions/binary/", this->charged_particle_list, this->target_particle_list);
    this->field_solver = read_voltage_inputs("../inputs/geometry.inp", this->scheme_type, *this->world);
    double plasma_freq = get_plasma_frequency(this->charged_particle_list[0].average_temperature, this->charged_particle_list[0].average_density);
    // Time step
    if (del_t_fraction / plasma_freq < del_t_temp) {
        this->del_t = del_t_fraction / plasma_freq;
        this->inv_plasma_freq_fraction = del_t_fraction;
    } else {
        this->del_t = del_t_temp;
        this->inv_plasma_freq_fraction = del_t_temp * plasma_freq;
    }
    if (this->del_t > this->simulation_time || this->simulation_time <= 0) {
        std::cout << "Issue with simulation time!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (this->averaging_time > 0 && this->del_t > this->averaging_time) {
        std::cout << "Issue with averaging time!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (this->number_diagnostics <= 0 || this->number_diagnostics > this->simulation_time/this->del_t) {
        std::cout << "Issue with number diagnostics!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "" << std::endl;
        std::cout << "Time setup:" << std::endl;
        std::cout << "------------------------------" << std::endl;
        std::cout << "Time step (s): " << this->del_t  << std::endl;
        std::cout << "Fraction of inverse plasma frequency: " << this->inv_plasma_freq_fraction << std::endl;
        std::cout << "Simulation time (s): " << this->simulation_time << std::endl;
        std::cout << "Averaging time (s): " << this->averaging_time << std::endl;
        std::cout << "Number diagnostics: " << this->number_diagnostics << std::endl;
        std::cout << "------------------------------" << std::endl;
        std::cout << " " << std::endl;
    }
    int save_str_len = this->save_file_folder.size();
    int dir_str_len = this->save_file_path.size();
    MPI_Bcast(&save_str_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dir_str_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_vars::mpi_rank != 0) {
        this->save_file_folder.resize(save_str_len);
        this->save_file_path.resize(dir_str_len);
    }
    MPI_Bcast(&this->save_file_folder[0], save_str_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&this->save_file_path[0], dir_str_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "" << std::endl;
        std::cout << "Simulation Directories:" << std::endl;
        std::cout << "------------------------------" << std::endl;
        std::cout << "Save folder name: " << this->save_file_folder << std::endl;
        std::cout << "Save folder path: " << this->save_file_path << std::endl;
        std::string restart;
        if (this->restarted_simulation) {
            restart = "Yes";
        } else {
            restart = "No";
        }
        std::cout << "Restarted simulation? " << restart << std::endl;
        std::cout << "------------------------------" << std::endl;
        std::cout << " " << std::endl;
    }
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++){
            this->charged_particle_list[part_num].sort_particle_diagnostics(thread_id, world->number_cells);
            this->charged_particle_list[part_num].gather_mpi();
        }
        #pragma omp barrier
        #pragma omp master
        {
            if (mpi_vars::mpi_rank == 0) {
                double sum = 0.0;
                for (int i = 0; i < this->charged_particle_list.size(); i++) {
                    sum += this->charged_particle_list[i].total_sum_v[0] * this->charged_particle_list[i].mass;
                }
                std::cout << "Total mv_x: " << sum << std::endl;
            }
        }
        #pragma omp barrier
        this->field_solver->deposit_charge_density(this->charged_particle_list, thread_id);
        #pragma omp master
        {
            this->field_solver->solve_potential(0.0, *this->world);
            this->field_solver->make_EField(*this->world);
        }
        #pragma omp barrier
        this->field_solver->push_particles(thread_id, this->del_t, this->charged_particle_list, *this->world);
        #pragma omp barrier
        for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++){
            this->charged_particle_list[part_num].sort_particle_diagnostics(thread_id, world->number_cells);
            this->charged_particle_list[part_num].gather_mpi();
        }
        #pragma omp barrier
        #pragma omp master
        {
            if (mpi_vars::mpi_rank == 0) {
                double sum = 0.0;
                for (int i = 0; i < this->charged_particle_list.size(); i++) {
                    sum += this->charged_particle_list[i].total_sum_v[0] * this->charged_particle_list[i].mass;
                }
                std::cout << "Total mv_x: " << sum << std::endl;
            }
        }
        #pragma omp barrier
    }
}
