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

    // Constructor declaration
    Domain(const std::string& filename);
    double get_x_from_xi(double xi);
    double get_xi_from_x(double xi);
    // Method declaration
    void readWorld();
};

#endif // DOMAIN_H
