#include "null_collision.h"
#include <regex>


void read_null_collision_inputs(const std::string& filename, const std::vector<Particle> &particle_list, const std::vector<Target_Particle> &target_particle_list){
    std::cout << " "<< std::endl;
    std::cout << "Reading null collision inputs "<< std::endl;
    std::cout << "---------------------------------------- "<< std::endl;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    int number_primary_particles = 0;
    std::vector<int> particle_primary_idx;
    std::vector<std::vector<int>> target_indices;
    std::vector<std::vector<std::vector<int>>> product_indices;
    std::vector<int> number_energy_points;
    std::vector<int> number_collisions_per_primary;
    std::vector<std::vector<int>> collision_type_per_primary;
    std::vector<std::vector<int>> product_idx_per_primary;
    std::vector<std::vector<std::vector<double>>> energy_arrays;
    std::vector<std::vector<std::vector<double>>> sigma_arrays;
    double E_threshold, temp_var, E_scaling, sigma_scaling;

    std::string line;
    std::getline(file, line);
    while (line.find("END") == std::string::npos) {
        std::istringstream iss(line);
        std::string coll_filename;
        iss >> coll_filename;
        coll_filename = "../../CollisionData/" + coll_filename;
        std::cout << "Open file: " << coll_filename << std::endl;
        std::ifstream coll_file(coll_filename);
        if (!coll_file) {
            std::cerr << "Error: Unable to open file" << coll_filename << std::endl;
            return;
        }
        std::getline(coll_file, line);
        while (line.find("END") == std::string::npos) {
            if (line.find("REACTION") != std::string::npos) {
                std::string reaction_string;
                std::getline(coll_file, line);
                iss.clear();
                iss.str(line);
                iss >> reaction_string;
                std::cout << reaction_string << std::endl;
                size_t arrow_pos = reaction_string.find("->");
                std::string reactant_string = reaction_string.substr(0,arrow_pos);
                std::string product_string = reaction_string.substr(arrow_pos+2);
                std::regex pattern(R"(\[([^\]]+)\])");
                std::smatch match;
                int primary_idx, target_idx; 
                int coll_idx;
                bool primary_exists = false;
                int i;
                // Get indices of the primary and target particle
                while (std::regex_search(reactant_string, match, pattern)){
                    
                    for (i=0;i<global_inputs::number_charged_particles;i++){
                        if (particle_list[i].name == match[1]){
                            primary_idx = i;     
                            break;
                        }
                    }
                    for (i=0;i<global_inputs::number_target_particles;i++){
                        if (target_particle_list[i].name == match[1]){
                            target_idx = i;
                            break;
                        }
                    }
                    reactant_string = match.suffix().str();
                }
            
                
           
                for (coll_idx=0;coll_idx<number_primary_particles;coll_idx++){
                    if (particle_primary_idx[coll_idx] == primary_idx){
                        primary_exists = true;
                        break;
                    }
                }
                if (!primary_exists){
                    number_primary_particles += 1;
                    particle_primary_idx.push_back(primary_idx);
                    target_indices.push_back({target_idx});
                    number_collisions_per_primary.push_back({1});
                }
                else {
                    target_indices[coll_idx].push_back(target_idx);
                    number_collisions_per_primary[coll_idx]++;
                }

                std::vector<int> local_product_indx;
                while (std::regex_search(product_string, match, pattern)){
                    for (i=0;i<global_inputs::number_charged_particles;i++){
                        if (particle_list[i].name == match[1]){
                            local_product_indx.push_back(i);
                        }
                    }
                    for (i=0;i<global_inputs::number_target_particles;i++){
                        if (target_particle_list[i].name == match[1]){
                            target_idx = i;
                        }
                    }
                    product_string = match.suffix().str();
                }
                local_product_indx.push_back(target_idx); // put target index as last
                
                std::getline(coll_file, line);
                iss.clear();
                iss.str(line);
                iss >> E_threshold >> temp_var;
                std::getline(coll_file, line);
                iss.clear();
                iss.str(line);
                iss >> E_scaling >> sigma_scaling;
                std::string coll_string;
                std::getline(coll_file, line);
                iss.clear();
                iss.str(line);
                iss >> coll_string;
                int collision_type = 0;
                if (coll_string == "ELASTIC") {
                    collision_type = 1;
                } else if (coll_string == "EXCITATION"){
                    collision_type = 2;
                } else if (coll_string == "IONIZATION"){
                    collision_type = 3;
                    // make sure primary particle is first, then followed by by-product
                    for (i=0;i<3;i++){
                        if (primary_idx == local_product_indx[i]){
                            // switch primary to front
                            local_product_indx[3] = local_product_indx[0];
                            local_product_indx[0] = primary_idx;
                            local_product_indx[i] = local_product_indx[3];
                            break;
                        }
                    }
                    for (i=1;i<3;i++){
                        if (particle_list[local_product_indx[i]].mass == Constants::electron_mass){
                            // switch electron to back
                            local_product_indx[3] = local_product_indx[2];
                            local_product_indx[2] = local_product_indx[i];
                            local_product_indx[i] = local_product_indx[3];
                            break;
                        }
                    }
                    local_product_indx.pop_back();
                } else if (coll_string == "CHARGEEXCHANGE"){
                    collision_type = 4;
                }
                if (!primary_exists){
                    collision_type_per_primary.push_back({collision_type});
                    product_indices.push_back({local_product_indx});
                }
                else {
                    collision_type_per_primary[coll_idx].push_back(collision_type);
                    product_indices[coll_idx].push_back(local_product_indx);
                }
                
                std::getline(coll_file, line);
                while (line.find("---------") == std::string::npos) {
                    std::getline(coll_file, line);
                }
                double energy_point, sigma_point;
                std::vector<double> local_energy_array;
                std::vector<double> local_sigma_array;
                std::getline(coll_file, line);
                while (line.find("---------") == std::string::npos) {
                    iss.clear();
                    iss.str(line);
                    iss >> energy_point >> sigma_point;
                    local_energy_array.push_back(energy_point * E_scaling);
                    local_sigma_array.push_back(sigma_point * sigma_scaling);
                    std::getline(coll_file, line);
                }
                // find which particle indices they match to
                
                if (!primary_exists){
                    energy_arrays.push_back({local_energy_array});
                    sigma_arrays.push_back({local_sigma_array});
                }
                else {
                    energy_arrays[coll_idx].push_back(local_energy_array);
                    sigma_arrays[coll_idx].push_back(local_sigma_array);
                }


            }
            std::getline(coll_file, line);
        }
        
        int k, l, u;
        std::cout << "number primary particles " << number_primary_particles << " is " << particle_primary_idx.size() << std::endl;
        for (k=0; k<number_primary_particles; k++){
            std::cout << "------------------------------ " << std::endl;
            std::cout << "Primary particle " << k << " is " << particle_list[particle_primary_idx[k]].name << std::endl;
            std::cout << "Amount of collisions is " << number_collisions_per_primary[k] << std::endl;
            for (l=0; l<collision_type_per_primary[k].size();l++){
                std::cout << "Collision type is " << collision_type_per_primary[k][l] << std::endl;
                std::cout << "Product indices are: " << std::endl;
                for (u=0; u < product_indices[k][l].size();u++){
                    std::cout << product_indices[k][l][u] << std::endl;  
                }
            }
            std::cout << "------------------------------ " << std::endl;
        }
        coll_file.close();
        std::getline(file, line);
    }
    file.close();
    // std::vector<Particle> particle_list;
    // int count_number_particles = 0;
    // std::string line;
    // // First check which particles exist and how many there are
    // while (std::getline(file, line)) {
    //     if (line.find("ELECTRONS") != std::string::npos){
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         std::istringstream iss(line);
    //         std::string name;
    //         std::cout << "Found electron"<< std::endl;
    //         uint32_t num_part_thread;
    //         size_t factor;
    //         iss >> name >> num_part_thread >> factor;
    //         std::cout << name << " " << num_part_thread << " " << factor << std::endl;
    //         std::getline(file, line);
    //         count_number_particles++;
    //     }

    //     if (line.find("IONS") != std::string::npos){
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         while (line.find("-------") == std::string::npos) {
    //             std::istringstream iss(line);
    //             std::string name;
    //             std::cout << "Found ions"<< std::endl;
    //             uint32_t num_part_thread;
    //             double mass_in, charge_in;
    //             size_t factor;
    //             iss >> name >> mass_in >> charge_in >> num_part_thread >> factor;
    //             std::cout << name << " " << " " << mass_in << " " << charge_in << " " << num_part_thread << " " << factor << std::endl;
    //             std::getline(file, line);
    //             count_number_particles++;
    //         }
    //     }
    // }
    
    // global_inputs::number_charged_particles = count_number_particles;
    // std::cout << global_inputs::number_charged_particles << std::endl;
    // particle_list.reserve(global_inputs::number_charged_particles);
    // // get electrons
    // file.clear();
    // file.seekg(0, std::ios::beg);
    // while (std::getline(file, line)) {
    //     if (line.find("ELECTRONS") != std::string::npos){
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         std::istringstream iss(line);
    //         std::string name;
    //         uint32_t num_part_thread;
    //         size_t factor;
    //         iss >> name >> num_part_thread >> factor;
    //         factor = factor * static_cast<size_t>(num_part_thread);
    //         particle_list.emplace_back(Constants::electron_mass, -Constants::elementary_charge, num_part_thread, factor, name);
    //         std::getline(file, line);
    //     }
    // }

    // file.clear();
    // file.seekg(0, std::ios::beg);
    // while (std::getline(file, line)) {
    //     if (line.find("IONS") != std::string::npos){
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         std::getline(file, line);
    //         while (line.find("-------") == std::string::npos) {
    //             std::istringstream iss(line);
    //             std::string name;
    //             uint32_t num_part_thread;
    //             double mass_in, charge_in;
    //             size_t factor;
    //             iss >> name >> mass_in >> charge_in >> num_part_thread >> factor;
    //             factor = factor * static_cast<size_t>(num_part_thread);
    //             mass_in = mass_in * Constants::mass_amu - charge_in * Constants::electron_mass;
    //             charge_in = charge_in * Constants::elementary_charge;
    //             particle_list.emplace_back(mass_in, charge_in, num_part_thread, factor, name);
    //             std::getline(file, line);
    //         }
    //     }
    // }
    // file.close();

    // int i;
    // double T_ave;
    // for (i = 0; i < global_inputs::number_charged_particles; i++){
    //     particle_list[i].initialize_weight(global_inputs::initial_density, world.L_domain);
    //     if (particle_list[i].mass == Constants::electron_mass) {
    //         T_ave = global_inputs::temp_electrons;
    //     }
    //     else {
    //         T_ave = global_inputs::temp_ions;
    //     }
    //     particle_list[i].initialize_rand_uniform(T_ave, world);
    //     std::cout << "Particle #: " << i << std::endl;
    //     std::cout << "Name: " << particle_list[i].name << std::endl;
    //     std::cout << "Mass: " << particle_list[i].mass << std::endl;
    //     std::cout << "Charge: " << particle_list[i].charge << std::endl;
    //     std::cout << "Weight: " << particle_list[i].weight << std::endl;
    //     std::cout << "Number of particles per thread: " << particle_list[i].number_particles[0] << std::endl;
    //     std::cout << "Allocated total amount of particles per thread " << particle_list[i].final_idx << std::endl;
    //     std::cout << "Averge KE: " << particle_list[i].get_KE_ave() << std::endl;
    //     std::cout << "----------------- " << std::endl;
    // }

    

}