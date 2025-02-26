#ifndef DOMAIN_H
#define DOMAIN_H

#include <iostream>
#include <vector>

class Domain {
public:
    // Grid properties
    std::vector<double> grid; 
    double del_X, L_domain;

    // Boundary conditions and threading indices
    std::vector<int> boundary_conditions;
    std::vector<std::vector<int>> thread_node_indx;
    int num_thread_node_indx;

    // Constructor declaration
    Domain(const std::string& filename);
    double get_x_from_xi(double xi) const;
    double get_xi_from_x(double xi) const;
    // Method declaration
    void readWorld();
};

#endif // DOMAIN_H
