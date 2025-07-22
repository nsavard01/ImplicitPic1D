#include "simulation/simulation.hpp"
#include "globals/plasma_functions.hpp"
#include <chrono>
#include <iomanip>
#include "globals/write_functions.hpp"



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

void simulation::initialize_diagnostic_files() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "Initializing diagnostic files in save directory... ";
        std::string folder_name = this->save_file_path + this->save_file_folder;
        bool dirExists = directoryExists(this->save_file_path);
        if (dirExists) {
            bool saveFileExists = directoryExists(folder_name);
            if (saveFileExists) {
                // Ask user for permission to overwrite existing directory
                std::cout << "Save directory " << folder_name << " already exists. Are you sure you want to continue (yes/no)? ";
                std::string userInput;
                std::cin >> userInput;
    
                if (userInput != "yes" && userInput != "Yes") {
                    std::cout << "You have decided to create a new directory for the save files." << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                removeDirectoryContents(folder_name);  // Remove old data
            } else {
                // Create the top-level directory
                if (!createDirectory(folder_name)) {
                    std::cerr << "Failed to create main directory: " << folder_name << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
    
            // Create the necessary directories
            if  (!createDirectory(folder_name + "/domain")  || !createDirectory(folder_name + "/phi") ||
                !createDirectory(folder_name + "/charged_particles") || !createDirectory(folder_name + "/target_particles")) {
                std::cerr << "Save directory not successfully created!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            for (int i = 0; i < this->charged_particle_list.size(); i++) {
                if (!createDirectory(folder_name + "/charged_particles/" + this->charged_particle_list[i].name)) {
                    std::cerr << "Save directory not successfully created!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                if (!createDirectory(folder_name + "/charged_particles/" + this->charged_particle_list[i].name + "/density")) {
                    std::cerr << "Save directory not successfully created!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                if (!createDirectory(folder_name + "/charged_particles/" + this->charged_particle_list[i].name + "/temperature")) {
                    std::cerr << "Save directory not successfully created!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                if (!createDirectory(folder_name + "/charged_particles/" + this->charged_particle_list[i].name + "/phase_space")) {
                    std::cerr << "Save directory not successfully created!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                if (this->null_collider_list[i].number_targets > 0) {
                    if (!createDirectory(folder_name + "/charged_particles/" + this->charged_particle_list[i].name + "/null_collision")) {
                        std::cerr << "Save directory not successfully created!" << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    for (int t_idx = 0; t_idx < this->null_collider_list[i].number_targets;t_idx++){
                        int idx = this->null_collider_list[i].target_idx[t_idx];
                        if (!createDirectory(folder_name + "/charged_particles/" + this->charged_particle_list[i].name + "/null_collision/" + this->target_particle_list[idx].name)) {
                            std::cerr << "Save directory not successfully created!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                    }
                }
            }

            for (int i = 0; i < this->target_particle_list.size(); i++) {
                if (!createDirectory(folder_name + "/target_particles/" + this->target_particle_list[i].name)) {
                    std::cerr << "Save directory not successfully created!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }


        
    
            // Copy input files (assuming the same mechanism for copying)
            std::string copyCommand = "cp -Tr ../inputs " + folder_name + "/inputs";
            int status = system(copyCommand.c_str());
            if (status != 0) {
                std::cerr << "Error copying input data deck" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        } else {
            std::cout << "Directory chosen to save data doesn't exist!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }


        // Initialize data files

        // Put Current date and time
        auto now = std::chrono::system_clock::now();
        std::time_t t = std::chrono::system_clock::to_time_t(now);
        
        // Convert to UTC time structure
        std::tm utc_time = *std::gmtime(&t);

        // Open file
        std::ofstream file(folder_name + "/date_time.dat");
        if (!file) {
            std::cerr << "Error opening file for date and time \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "UTC_Date UTC_Time\n";

        // Write date and time in YYYY MM DD HH MM SS format
        file << std::put_time(&utc_time, "%Y%m%d %H%M%S") << "\n";

        file.close();

        file.open(folder_name + "/initial_condition.dat");

        // Write header (optional)
        file << "number MPI, number threads, scheme_type, Final Expected Time(s), Delta t(s), numDiag \n";
        file << std::scientific << std::setprecision(8);
        file << mpi_vars::mpi_size << "\t"
            << omp_get_max_threads() << "\t"
            << this->scheme_type << "\t"
            << this->simulation_time << "\t"
            << this->del_t << "\t"
            << this->number_diagnostics << "\t"
            << this->averaging_time
            <<"\n";

        file.close();

        file.open(folder_name + "/simulation_timing_data.dat");

        // Write header (optional)
        file << "Total time (s), field time (s), particle_time (s), null collision time (s), artificial collision time (s) \n";

        file.close();

        file.open(folder_name + "/global_diagnostic_data.dat");

        // Write header (optional)
        file << "Time (s), total steps, total momentum x (kg /s /m), total momentum y (kg /s /m), total momentum z (kg /s /m), total energy (J/m^2) \n";

        file.close();

        this->world->write_domain(folder_name);
        this->field_solver->initialize_diagnostic_files(folder_name);

        for (int part_num = 0; part_num < this->target_particle_list.size(); part_num++) {
            this->target_particle_list[part_num].initialize_diagnostic_files(folder_name);
        }
        for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++) {
            this->charged_particle_list[part_num].initialize_diagnostic_files(folder_name);
            this->null_collider_list[part_num].initialize_diagnostic_files(folder_name, this->charged_particle_list, this->target_particle_list);
        }

        std::cout << "Done" << std::endl;

    }
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
            std::cout << "Scheme type: 1 (EC-PIC)" << std::endl;
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
    for (int i = 0; i < charged_particle_list.size(); i++){
        this->null_collider_list[i].set_null_frequency(this->charged_particle_list, this->target_particle_list);
        this->null_collider_list[i].print_out(this->charged_particle_list, this->target_particle_list);
    }
    this->particle_operator_list = read_particle_operators("../inputs/particle_operations/", this->charged_particle_list, *this->world);
    for (int i = 0; i < particle_operator_list.size(); i++){
        particle_operator_list[i]->print_out();
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
    this->field_solver = read_voltage_inputs("../inputs/geometry.inp", this->scheme_type, *this->world);
    double plasma_freq = get_plasma_frequency(this->charged_particle_list[0].average_density);
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
    if (this->number_diagnostics <= 0 || this->number_diagnostics > std::ceil(this->simulation_time/this->del_t) + 1) {
        std::cout << "Issue with number diagnostics! " << std::endl;
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
    // Setup directories
    this->initialize_diagnostic_files();
    this->save_file_folder = this->save_file_path + this->save_file_folder;
    this->simulation_start_time = 0.0;
    this->diag_time_division = (this->simulation_time - this->simulation_start_time)/(this->number_diagnostics-1);
    this->current_diag_step = 0;
    this->current_time = 0.0;
    this->next_diag_time = this->current_time;
    this->last_diag_time = this->current_time;
    this->elapsed_time = 0.0;
    this->current_step = 0;
    this->particle_time = 0;
    this->field_time = 0;
    this->null_collision_time = 0;
    this->art_collision_time = 0;

    // Initialize potential and diagnostics
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "" << std::endl;
        std::cout << "Initializing potential and let us run!" << std::endl;
        std::cout << "------------" << std::endl;
     }
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        this->field_solver->deposit_charge_density(*this->world, this->charged_particle_list, thread_id);
        #pragma omp barrier
        #pragma omp master
        {
            this->field_solver->solve_potential(this->current_time, *this->world);
            this->field_solver->make_EField(*this->world);
        }
        #pragma omp barrier
        this->diagnostics(thread_id);
        #pragma omp barrier
        this->reset_diagnostics(thread_id);
        
    }
    

    

    

}

void simulation::diagnostics(int thread_id) {

    #pragma omp barrier
    // write initial diagnostics
    
    
    for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++){
        this->charged_particle_list[part_num].get_particle_diagnostics(thread_id, this->world->number_cells, 1);
    }
    // #pragma omp barrier
    // this->field_solver->deposit_density(this->charged_particle_list, thread_id);
    #pragma omp barrier
    #pragma omp master
    {
        double total_momentum[3];
        total_momentum[0] = total_momentum[1] = total_momentum[2] = 0.0;
        this->field_solver->get_diagnostics(*this->world, this->charged_particle_list);
        this->field_solver->write_diagnostics(this->save_file_folder, this->current_diag_step);
        double total_Energy = this->field_solver->total_field_energy;
        for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++){
            this->charged_particle_list[part_num].gather_mpi();
            this->charged_particle_list[part_num].write_diagnostics(this->save_file_folder, this->current_diag_step);
            total_momentum[0] += this->charged_particle_list[part_num].weight * this->charged_particle_list[part_num].mass * this->charged_particle_list[part_num].total_sum_v[0];
            total_momentum[1] += this->charged_particle_list[part_num].weight * this->charged_particle_list[part_num].mass * this->charged_particle_list[part_num].total_sum_v[1];
            total_momentum[2] += this->charged_particle_list[part_num].weight * this->charged_particle_list[part_num].mass * this->charged_particle_list[part_num].total_sum_v[2];
            double v_sqr = 0.0;
            for (int i = 0; i<3; i++) {
                v_sqr += this->charged_particle_list[part_num].total_sum_v_square[i];
            }
            total_Energy += 0.5 * this->charged_particle_list[part_num].weight * this->charged_particle_list[part_num].mass * v_sqr;
            this->null_collider_list[part_num].gather_mpi();
            this->null_collider_list[part_num].write_diagnostics(this->save_file_folder, this->charged_particle_list, this->target_particle_list);
        }
        this->field_solver->write_particle_densities(this->save_file_folder, "density_" + std::to_string(this->current_diag_step) + ".dat", this->charged_particle_list, *this->world);
        for (int t_idx = 0; t_idx < this->target_particle_list.size(); t_idx++) {
            this->target_particle_list[t_idx].write_diagnostics(this->save_file_folder, this->current_diag_step);
        }
        if (mpi_vars::mpi_rank == 0) {
            std::ofstream file(this->save_file_folder + "/global_diagnostic_data.dat", std::ios::app);

            file << std::scientific << std::setprecision(8);
            file << this->current_time << "\t"
            << this->current_step << "\t"
            << total_momentum[0] << "\t"
            << total_momentum[1] << "\t"
            << total_momentum[2] << "\t"
            << total_Energy << "\n";

            file.close();

            file.open(this->save_file_folder + "/simulation_timing_data.dat", std::ios::app);

            file << std::scientific << std::setprecision(8);
            file << this->elapsed_time << "\t"
            << this->field_time << "\t"
            << this->particle_time << "\t"
            << this->null_collision_time << "\t"
            << this->art_collision_time << "\n";

            file.close();

            std::cout << "" << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "Diagnostics # " << std::to_string(this->current_diag_step+1) << "/" << std::to_string(this->number_diagnostics) << std::endl;
            std::cout << "Current simulation time (s): " << this->current_time << std::endl;
            for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++){
                std::cout << "Particle " << this->charged_particle_list[part_num].name << ", total #: " << this->charged_particle_list[part_num].total_number_particles << " , n_ave: " 
                << this->charged_particle_list[part_num].average_density << " (1/m^3) , T_ave: " << this->charged_particle_list[part_num].average_temperature << " (eV) " << std::endl;
                if (this->null_collider_list[part_num].number_targets > 0) {
                    std::cout << "Null collisions for particle " << this->charged_particle_list[part_num].name << ": " << std::endl;
                    for (int t_idx = 0; t_idx < this->null_collider_list[part_num].number_targets; t_idx++) {
                        int idx = this->null_collider_list[part_num].target_idx[t_idx];
                        double freq = 0.0;
                        for (int coll_idx = 0; coll_idx < this->null_collider_list[part_num].number_collisions_per_target[t_idx]; coll_idx++) {
                            freq += double(this->null_collider_list[part_num].total_amount_collisions[t_idx][coll_idx]);
                        }
                        freq = freq / double(this->null_collider_list[part_num].total_amount_collidable_particles) / this->del_t;
                        std::cout << "Target particle: " << this->target_particle_list[idx].name << ", collision frequency (Hz): " 
                        << freq << std::endl;
                    }
                }
            }
            std::cout << "Total momentum (kg / m/s) is x: " << total_momentum[0] << " y: " << total_momentum[1] << " z: " << total_momentum[2] << std::endl;
            std::cout << "Total Energy (J/m^2) is : " << total_Energy << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "" << std::endl;
        }
    }

    #pragma omp barrier
}

void simulation::reset_diagnostics(int thread_id) {
    for (int i = 0; i<this->charged_particle_list.size(); i++){
        this->charged_particle_list[i].reset_diagnostics(thread_id);
        this->null_collider_list[i].order_collisions();
        this->null_collider_list[i].reset_diagnostics(thread_id);
    }
    #pragma omp master
    {   
        this->current_diag_step++;
        this->last_diag_time = this->next_diag_time;
        this->next_diag_time = this->last_diag_time + this->diag_time_division;
        this->diag_step_diff = 0;
    }
}

void simulation::run() {

    double start_time_total = MPI_Wtime();
    double timer_1, timer_2;
    int number_charged_particles = this->charged_particle_list.size();

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        while (this->current_time < this->simulation_time) {
            this->field_solver->integrate_time_step(thread_id, this->del_t, this->current_time, *this->world, this->charged_particle_list);
            #pragma omp barrier
            #pragma omp master
            {
                this->particle_time += this->field_solver->particle_timer;
                this->field_time += this->field_solver->potential_timer;
                
                timer_1 = MPI_Wtime();  
            }
            for (int part_num = 0; part_num < number_charged_particles; part_num++){
                this->null_collider_list[part_num].generate_null_collisions(thread_id, this->charged_particle_list, this->target_particle_list, this->del_t);
            }
            #pragma omp barrier
            #pragma omp master
            {   
                timer_2 = MPI_Wtime();
               
                this->null_collision_time += (timer_2 - timer_1);
                this->current_time += this->del_t;
                this->current_step++;
                this->diag_step_diff++;
                // if (mpi_vars::mpi_rank == 0) {
                //     if (this->current_step > 10) {
                //         MPI_Abort(MPI_COMM_WORLD, 1);
                //     }
                //     std::cout << "Current time (s): " << this->current_time << ", Current step: " << this->current_step << std::endl;
                // }
            }
            #pragma omp barrier
            if (this->current_time >= this->next_diag_time) {
                #pragma omp master
                {
                    
                    this->elapsed_time = MPI_Wtime() - start_time_total;
                }
                this->diagnostics(thread_id);
                #pragma omp barrier
                this->reset_diagnostics(thread_id);
            }
            #pragma omp barrier
        }
    }
    double end_time_total = MPI_Wtime();
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "Simulation took " << end_time_total - start_time_total << " seconds" << std::endl;
    }

}

void simulation::averaging() {

    if (this->averaging_time > 0) {
        if (mpi_vars::mpi_rank == 0) {
            std::cout << "" << std::endl;
            std::cout << "Averaging diagnostics over time: " << this->averaging_time << " seconds" << std::endl;
            std::cout << "----------------------------------- " << std::endl;
            std::cout << "" << std::endl;


        }

        
        double end_simulation_time = this->current_time + this->averaging_time;
        int number_charged_particles = this->charged_particle_list.size();
        double start_sim_time = this->current_time;
        this->current_step = 0;
        double res = 1.0;
        if (this->field_solver->RF_rad_frequency > 0.0) {
            double RF_period = 2.0 * M_PI / this->field_solver->RF_rad_frequency;
            this->diag_time_division = RF_period; // 1 RF periods
        } else {
            if (this->charged_particle_list.size() > 0 && this->charged_particle_list[0].mass == constants::electron_mass) {
                double n_e = this->charged_particle_list[0].average_density/this->world->length_domain;
                double plasma_freq = get_plasma_frequency(n_e);
                this->diag_time_division = 50.0 / plasma_freq; // 50 plasma periods
            } else {
                this->diag_time_division = 100.0 * this->del_t; // Default value 100 * del_t
            }
        }
        this->next_diag_time = this->current_time + this->diag_time_division;
        std::vector<double> average_phi = this->field_solver->phi;
        std::vector<double> average_phi_check = this->field_solver->phi;

        // Initialize vectors for distribution functions
        double start_time = MPI_Wtime();
        
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++){
                this->charged_particle_list[part_num].get_particle_diagnostics(thread_id, this->world->number_cells, 1);
            }
            #pragma omp barrier
            while (this->current_time < end_simulation_time && res > 1e-6) {
                this->field_solver->integrate_time_step(thread_id, this->del_t, this->current_time, *this->world, this->charged_particle_list);
                #pragma omp barrier
                #pragma omp for
                for (int i = 0; i < this->world->number_nodes; i++) {
                    average_phi[i] += this->field_solver->phi[i];
                }
                #pragma omp barrier
                for (int part_num = 0; part_num < this->charged_particle_list.size(); part_num++){
                    this->charged_particle_list[part_num].get_particle_diagnostics(thread_id, this->world->number_cells, 1);
                }
                #pragma omp barrier
                for (int part_num = 0; part_num < number_charged_particles; part_num++){
                    this->null_collider_list[part_num].generate_null_collisions(thread_id, this->charged_particle_list, this->target_particle_list, this->del_t);
                }
                #pragma omp barrier
                #pragma omp master
                {   
                    this->current_time += this->del_t;
                    this->current_step++;
                }
                #pragma omp barrier
                if (this->current_time >= this->next_diag_time) {
                    // average phi and check if can break loop early
                    #pragma omp master
                    {
                        res = 0.0;
                        double curr_average;
                        double integ_average = 0.0;
                        double diff;
                        for (int i = 0; i < this->world->number_nodes; i++) {
                            // compute current running average
                            curr_average = average_phi[i] / double(this->current_step + 1);
                            // compute difference with previous average
                            integ_average += curr_average;
                            diff = average_phi_check[i] - curr_average;
                            // compute difference percentage squared
                            res += diff*diff;
                            // update average_phi_check
                            average_phi_check[i] = curr_average;
                        }
                        res = std::sqrt(res / double(this->world->number_nodes));
                        res = res/integ_average;
                        this->next_diag_time = this->current_time + this->diag_time_division;
                    }
                }
                #pragma omp barrier
            }
            #pragma omp barrier
        }
        double end_time = MPI_Wtime();
        for (int part_num = 0; part_num < number_charged_particles; part_num++){
            this->charged_particle_list[part_num].gather_mpi();
        }
        if (mpi_vars::mpi_rank == 0) {
            write_vector_to_binary_file(average_phi_check, this->world->number_nodes, this->save_file_folder + "/phi/potential_average.dat", 0);
            std::cout << "out here" << std::endl;
            for (int part_num = 0; part_num < number_charged_particles; part_num++){
                this->charged_particle_list[part_num].write_diagnostics_average(this->save_file_folder);
                for (int i = 0; i < this->world->number_nodes; i++) {
                    this->charged_particle_list[part_num].density[i] /= double(this->current_step + 1);
                }
            }
            std::cout << "out here again" << std::endl;
            this->field_solver->write_particle_densities(this->save_file_folder, "density_average.dat", this->charged_particle_list, *this->world);
            std::cout << "Averaging finished and took " << end_time - start_time <<  " seconds" << std::endl;
            std::cout << "Ended over simulation time of " << this->current_time - start_sim_time << std::endl;
            std::cout << "Final res is: " << res << std::endl;
        }
        

    } else {
        std::cout << "" << std::endl;
        std::cout << "No averaging!" << std::endl;
        std::cout << "----------------------------------- " << std::endl;
        std::cout << "" << std::endl;
    }
    

}
