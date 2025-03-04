#include "simulation.h"
#include <iostream>
#include "basic_tools.h"
#include <cstdlib>
#include <mpi.h>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <chrono>





// Function to remove a directory and its contents
void removeDirectoryContents(const std::string& dirName) {
    std::string command = "rm -r " + dirName + "/*";
    int status = system(command.c_str());
    if (status != 0) {
        std::cerr << "Error removing directory contents: " << dirName << std::endl;
        exit(1);  // Exiting since the function isn't able to clean up properly
    }
}

Simulation::Simulation(const std::string dirName){
    if (Constants::mpi_rank == 0) {
        bool dirExists = directoryExists(dirName);
        
        if (dirExists) {
            // Ask user for permission to overwrite existing directory
            std::cout << "Save directory " << dirName << " already exists. Are you sure you want to continue (yes/no)? ";
            std::string userInput;
            std::cin >> userInput;

            if (userInput != "yes" && userInput != "Yes") {
                std::cout << "You have decided to create a new directory for the save files." << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            removeDirectoryContents(dirName);  // Remove old data
        } else {
            // Create the top-level directory
            if (!createDirectory(dirName)) {
                std::cerr << "Failed to create main directory: " << dirName << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Create the necessary directories
        if (!createDirectory(dirName + "/Density") ||
            !createDirectory(dirName + "/PhaseSpace") ||
            !createDirectory(dirName + "/Phi") ||
            !createDirectory(dirName + "/Temperature")) {
            std::cerr << "Save directory not successfully created!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Create BinaryCollisions directory if collNumber > 0
        if (global_inputs::number_binary_collisions > 0) {
            if (!createDirectory(dirName + "/BinaryCollisions")) {
                std::cerr << "Save directory not successfully created!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Copy input files (assuming the same mechanism for copying)
        std::string copyCommand = "cp -Tr ../InputData " + dirName + "/InputData";
        int status = system(copyCommand.c_str());
        if (status != 0) {
            std::cerr << "Error copying input data deck" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    this->directory_name = dirName;

}

void Simulation::initialize_data_files(Potential_Solver& solver, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list,
    std::vector<Null_Collision>& binary_collision_list, const Domain& world){
    
    if (Constants::mpi_rank == 0) {
        // Initialize files for data dumps

        // Put Current date and time
        auto now = std::chrono::system_clock::now();
        std::time_t t = std::chrono::system_clock::to_time_t(now);
        
        // Convert to UTC time structure
        std::tm utc_time = *std::gmtime(&t);

        // Open file
        std::ofstream file(this->directory_name + "/DateTime.dat");
        if (!file) {
            std::cerr << "Error opening file for DateTime \n";
            return;
        }

        // Write header (optional)
        file << "UTC_Date UTC_Time\n";

        // Write date and time in YYYY MM DD HH MM SS format
        file << std::put_time(&utc_time, "%Y%m%d %H%M%S") << "\n";

        file.close();

        file.open(this->directory_name + "/InitialConditions.dat");

        // Write header (optional)
        file << "Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, delX, n_ave, T_e, T_i, numDiag, NumChargedPart, numThread, RF_frequency, RF_half_amplitude \n";
        file << std::scientific << std::setprecision(8);
        // Write date and time in YYYY MM DD HH MM SS format
        file << global_inputs::number_nodes << "\t"
            << global_inputs::simulation_time << "\t"
            << global_inputs::time_step << "\t"
            << global_inputs::time_step << "\t"
            << world.del_X << "\t"
            << global_inputs::initial_density << "\t"
            << global_inputs::temp_ions << "\t"
            << global_inputs::temp_electrons << "\t"
            << global_inputs::number_diagnostic_steps << "\t"
            << global_inputs::number_charged_particles << "\t"
            << global_inputs::number_omp_threads << "\t"
            << solver.RF_rad_frequency << "\t"
            << solver.RF_half_amplitude
            <<"\n";

        file.close();

        file.open(this->directory_name + "/ParticleProperties.dat");

        // Write header (optional)
        file << "Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2), maxIdx \n";
        file << std::scientific << std::setprecision(8);
        
        for (int i = 0;i<global_inputs::number_charged_particles;i++){
            file << particle_list[i].name << "\t"
                << particle_list[i].mass << "\t"
                << particle_list[i].charge << "\t"
                << particle_list[i].weight << "\t"
                << particle_list[i].final_idx
                <<"\n";
        }
        

        file.close();

        file.open(this->directory_name + "/TargetParticleProperties.dat");

        // Write header (optional)
        file << "Particle Symbol, Particle Mass (kg), Particle Temp (K), Particle Density (m^-3) \n";
        file << std::scientific << std::setprecision(8);
        
        for (int i = 0;i<global_inputs::number_target_particles;i++){
            file << target_particle_list[i].name << "\t"
                << target_particle_list[i].mass << "\t"
                << target_particle_list[i].temperature << "\t"
                << target_particle_list[i].density
                <<"\n";
        }
        

        file.close();

        file.open(this->directory_name + "/SimulationTimeData.dat");

        // Write header (optional)
        file << "Elapsed Times(s), Potential Time (s), Mover Time (s), Collision Time (s), Total Steps \n";

        file.close();

        file.open(this->directory_name + "/SimulationTimeData.dat");

        // Write header (optional)
        file << "Elapsed Times(s), Potential Time (s), Mover Time (s), Collision Time (s), Total Steps \n";

        file.close();

        file.open(this->directory_name + "/GlobalDiagnosticData.dat");

        // Write header (optional)
        file << "Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), TotalMomentum(kg/m/s), TotalEnergy(J/m^2) \n";

        file.close();

        file.open(this->directory_name + "/SimulationFinalData.dat");

        // Write header (optional)
        file << "Elapsed Times(s), Potential Time (s), Mover Time (s), Collision Time (s), Total Steps \n";

        file.close();

        
        for (int i = 0;i<global_inputs::number_charged_particles;i++){
            particle_list[i].initialize_diagnostic_file(this->directory_name);
        }

        for (int i=0; i< global_inputs::number_binary_collisions; i++){
            binary_collision_list[i].initialize_data_files(this->directory_name, particle_list, target_particle_list);
        }
        

    }


}

void Simulation::reset_diag(std::vector<Particle>& particle_list, std::vector<Null_Collision>& binary_collision_list) {
    for (int i=0;i<global_inputs::number_charged_particles;i++){
        particle_list[i].accum_wall_loss[0] = 0;
        particle_list[i].accum_wall_loss[1] = 0;
        particle_list[i].accum_energy_loss[0] = 0.0;
        particle_list[i].accum_energy_loss[1] = 0.0;
        particle_list[i].accum_momentum_loss[0] = 0.0;
        particle_list[i].accum_momentum_loss[1] = 0.0;
    }

    for (int i=0;i<global_inputs::number_binary_collisions;i++){
        binary_collision_list[i].total_amount_collidable_particles = 0;
        for (int j = 0;j<binary_collision_list[i].number_collisions;j++){     
            binary_collision_list[i].total_incident_energy[j] = 0.0;
            binary_collision_list[i].total_energy_loss[j] = 0.0;
            binary_collision_list[i].total_amount_collisions[j] = 0.0;
        }
    }

    this->charge_loss = 0;
    this->energy_loss = 0;
    this->collision_energy_loss = 0;
    this->momentum_loss = 0;
    this->momentum_total[0] = this->momentum_total[1] = this->momentum_total[2] = 0;

    this->diag_step_diff = 0;
}

void Simulation::diag_write(Potential_Solver& solver, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list,
    std::vector<Null_Collision>& binary_collision_list, const Domain& world) {
    
    double time_diff = this->diag_step_diff * global_inputs::time_step;
    solver.write_phi(this->directory_name, this->current_diag_step, false);
    double total_energy_system = solver.get_total_PE(world);
    for (int  i=0;i<global_inputs::number_charged_particles;i++){
        particle_list[i].gather_mpi();
        particle_list[i].load_density(true);
        particle_list[i].write_density(this->directory_name, world, this->current_diag_step, false);
        particle_list[i].write_cell_temperature(this->directory_name, this->current_diag_step);
        particle_list[i].diag_write(this->directory_name, time_diff, this->current_time);
        particle_list[i].write_phase_space(this->directory_name, this->current_diag_step);
        total_energy_system += particle_list[i].get_KE_total();

        this->charge_loss += (particle_list[i].accum_wall_loss[0] + particle_list[i].accum_wall_loss[1]) * particle_list[i].q_times_wp; 
        this->energy_loss += (particle_list[i].accum_energy_loss[0] + particle_list[i].accum_energy_loss[1]) * particle_list[i].mass * 0.5 * particle_list[i].weight;
        this->momentum_loss += (particle_list[i].accum_momentum_loss[0] + particle_list[i].accum_momentum_loss[1]) * particle_list[i].mass * particle_list[i].weight;
        this->momentum_loss += particle_list[i].get_momentum_total() * particle_list[i].weight;
    }
    for (int  i=0;i<global_inputs::number_binary_collisions;i++){
        binary_collision_list[i].gather_mpi();
        binary_collision_list[i].diag_write(this->directory_name, particle_list, target_particle_list, time_diff);
        for (int coll_num = 0;coll_num < binary_collision_list[i].number_collisions;coll_num++){
            this->collision_energy_loss += 0.5 * binary_collision_list[i].total_energy_loss[coll_num] * particle_list[binary_collision_list[i].primary_idx].weight;
        }
    }

    if (Constants::mpi_rank == 0){
        std::ofstream file(this->directory_name + "/GlobalDiagnosticData.dat", std::ios::app); 
        file << std::scientific << std::setprecision(8);
        file << current_time << "\t"
            << this->collision_energy_loss/time_diff << "\t"
            << this->charge_loss/time_diff << "\t"
            << this->energy_loss/time_diff << "\t"
            << this->momentum_loss << "\t"
            << total_energy_system
            <<"\n";
        file.close();

        this->end_time_total = MPI_Wtime();

        file.open(this->directory_name + "/SimulationTimeData.dat", std::ios::app); 
        file << std::scientific << std::setprecision(8);
        file << this->end_time_total - this->start_time_total << "\t"
            << this->tot_potential_time << "\t"
            << this->tot_particle_time << "\t"
            << this->tot_collision_time << "\t"
            << this->current_step
            <<"\n";
        file.close();

    }
    
}

void Simulation::run(Potential_Solver& solver, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list,
    std::vector<Null_Collision>& binary_collision_list, const Domain& world) {


    // Take initial diagnostics
    world.write_domain(this->directory_name);
    for (int  i=0;i<global_inputs::number_charged_particles;i++){
        particle_list[i].load_density(true);
        particle_list[i].write_density(this->directory_name, world, 0, false);
        particle_list[i].write_cell_temperature(this->directory_name, 0);
    }
    solver.write_phi(this->directory_name, 0, false);

    
    this->current_time = global_inputs::start_simulation_time;
    this->diag_time_division = (global_inputs::simulation_time - this->current_time)/global_inputs::number_diagnostic_steps;
    this->diag_time = this->current_time + this->diag_time_division;
    this->current_diag_step = 1;
    this->elapsed_time = 0;
    this->diag_step_diff = 0;

    this->start_time_total = MPI_Wtime();
    double start_time;
    double end_time;
    this->tot_particle_time = 0.0;
    this->tot_potential_time = 0.0;
    this->tot_collision_time = 0.0;
    this->current_step = 0;
   
    while (this->current_time < global_inputs::simulation_time) { 
        this->current_time += global_inputs::time_step;
        start_time = MPI_Wtime();
        solver.move_particles(particle_list, world, global_inputs::time_step);
        solver.deposit_rho(particle_list, world);
        end_time = MPI_Wtime();
        this->tot_particle_time += (end_time - start_time);

        start_time = MPI_Wtime();
        solver.solve_potential_tridiag(world, this->current_time);
        solver.make_EField(world);
        end_time = MPI_Wtime();
        this->tot_potential_time += (end_time - start_time);

        start_time = MPI_Wtime();
        for (int coll_num = 0;coll_num < global_inputs::number_binary_collisions;coll_num++){
            binary_collision_list[coll_num].generate_null_collisions(particle_list, target_particle_list, global_inputs::time_step);
        }
        end_time = MPI_Wtime();
        this->tot_collision_time += (end_time - start_time);

        this->diag_step_diff++;
        this->current_step++;

        if (this->current_time >= this->diag_time) {
            if (Constants::mpi_rank == 0) {
                std::cout << "Simulation is " << this->current_time * 100 /global_inputs::simulation_time << " % done" << std::endl;
            }
            this->diag_write(solver, particle_list, target_particle_list, binary_collision_list, world);
            this->reset_diag(particle_list, binary_collision_list);
            this->diag_time += this->diag_time_division;
            this->current_diag_step++;
        }
    }
    this->end_time_total = MPI_Wtime();
    if (Constants::mpi_rank == 0) {
        std::ofstream file(this->directory_name + "/SimulationFinalData.dat", std::ios::app);
        file << std::scientific << std::setprecision(8);
        file << this->end_time_total - this->start_time_total << "\t"
            << this->tot_potential_time << "\t"
            << this->tot_particle_time << "\t"
            << this->tot_collision_time << "\t"
            << this->current_step
            <<"\n";
        file.close();
        std::cout << "Total particle time " << this->tot_particle_time << " seconds" << std::endl;
        std::cout << "Total potential time " << this->tot_potential_time << " seconds" << std::endl;
        std::cout << "Total collision time " << this->tot_collision_time << " seconds" << std::endl;
        std::cout << "Total time is " << this->end_time_total - this->start_time_total << " seconds" << std::endl;
    }
    
}




