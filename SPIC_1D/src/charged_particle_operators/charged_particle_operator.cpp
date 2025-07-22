
#include "charged_particle_operators/charged_particle_operator.hpp"
#include "charged_particle_operators/charged_particle_wall_injector.hpp"
#include <iomanip>
#include <dirent.h>
#include <fstream>
#include <sstream>


std::vector<std::unique_ptr<charged_particle_operator>> read_particle_operators(const std::string& directory_path, const std::vector<charged_particle>& particle_list, const domain& world) {
    std::vector<std::unique_ptr<charged_particle_operator>> output;
    bool wall_injection_bool = false;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "<< std::endl;
        std::cout << "Reading charged particle operations "<< std::endl;
        std::cout << "---------------------------------------- "<< std::endl;
    }
    for (int rank_num = 0; rank_num < mpi_vars::mpi_size; rank_num++){
        if (mpi_vars::mpi_rank == rank_num) {
            // Open the directory
            DIR* dir = opendir(directory_path.c_str());
            if (!dir) {
                perror("opendir");
                exit(EXIT_FAILURE);
            }

            struct dirent* entry;
            while ((entry = readdir(dir)) != nullptr) {
                if (entry->d_type == DT_REG) {  // regular file
                    std::string filename = directory_path + entry->d_name;
                    if (!filename.compare(directory_path + "wall_injection.inp")) {
                        wall_injection_bool = true;
                        std::vector<double> current_density, v_x_array;
                        std::vector<int> direction, wall_node_location, particle_indx;
                        std::string filename = directory_path + entry->d_name;
                        std::string line;
                        std::ifstream file(filename);
                        if (!file) {
                            std::cerr << "Error: Unable to open file " << filename << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        std::getline(file, line);
                        std::istringstream iss(line);
                        while (line.find("END") == std::string::npos) {
                            if (line.find("----") != std::string::npos) {
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                std::string name;
                                iss >> name;
                                int indx = -1;
                                for (int i = 0; i < particle_list.size(); i++) {
                                    if (particle_list[i].name == name) {
                                        indx = i;
                                        break;
                                    }
                                }
                                if (indx < 0) {
                                    std::cout << "particle " << name << " does not exist! " << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                particle_indx.push_back(indx);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                int node;
                                iss >> node;
                                if (node != 0 && node != world.number_nodes) {
                                    std::cout << "WARNING: Node for particle injection not on boundary!" << std::endl;
                                } else if (node == 0) {
                                    if ((world.left_boundary_condition != 1) && (world.left_boundary_condition != 4)){
                                        std::cout << "WARNING: Leftmost node for particle injection not on metallic boundary!" << std::endl;
                                    }
                                } else if (node == world.number_nodes) {
                                    if ((world.right_boundary_condition != 1) && (world.right_boundary_condition != 4)){
                                        std::cout << "WARNING: Rightmost node for particle injection not on metallic boundary!" << std::endl;
                                    }
                                }
                                wall_node_location.push_back(node);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double J;
                                iss >> J;
                                J = std::abs(J);
                                current_density.push_back(J);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double v_x;
                                iss >> v_x;
                                v_x_array.push_back(v_x);
                                iss.clear();
                                std::getline(file, line);
                            }
                            
                            std::getline(file, line);
                        }
                        file.close();
                        output.push_back(std::make_unique<charged_particle_wall_injector>(particle_indx, current_density, v_x_array, wall_node_location));
                    }
                }
            }   
        
            closedir(dir);
        
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return output;

}
