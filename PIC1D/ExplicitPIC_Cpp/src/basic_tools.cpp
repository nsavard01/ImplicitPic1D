#include "basic_tools.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>
#include "Constants.h"
#include <cmath>

template <typename T>
void write_vector_to_binary_file(const std::vector<T>& vec, const std::string& filename, size_t paddingBytes, char paddingValue, bool append) {
    // Calculate total size (data + padding)
    size_t size_of_T = sizeof(T);
    size_t total_size = 2* paddingBytes + vec.size() * size_of_T;
    

    // Create a byte buffer and pre-allocate required space
    std::vector<char> byte_stream(total_size);

    // Add padding if required
    if (paddingBytes > 0) {
        std::fill(byte_stream.begin(), byte_stream.begin() + paddingBytes, paddingValue);  // Fill padding at the beginning
    }

    // Copy the vector data into the byte buffer (after padding)
    for (size_t i = 0; i < vec.size(); ++i) {
        std::memcpy(byte_stream.data() + paddingBytes + i * size_of_T, &vec[i], size_of_T);
    }

    if (paddingBytes > 0) {
        std::fill(byte_stream.begin() + paddingBytes + vec.size() * size_of_T, byte_stream.end(), paddingValue);  // Fill padding at the end
    }

    std::ios_base::openmode mode = std::ios::binary;
    if (append){
        mode |= std::ios::app;
    }

    // Write the byte buffer to the binary file
    std::ofstream outFile(filename, mode);
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    outFile.write(byte_stream.data(), byte_stream.size());  // Write raw bytes to file
    outFile.close();

}

// Function to create a directory
bool createDirectory(const std::string& dirName) {
    if (mkdir(dirName.c_str(), 0777) == -1) {
        std::cerr << "Error creating directory: " << dirName << std::endl;
        return false;
    }
    return true;
}

bool directoryExists(const std::string& dirName) {
    DIR* dir = opendir(dirName.c_str());
    if (dir) {
        closedir(dir);
        return true;
    }
    return false;
}

double plasma_frequency(double n_e) {
    return Constants::elementary_charge * std::sqrt(n_e /Constants::electron_mass / Constants::epsilon_0);
}

double debye_length(double n_e, double T_e) {
    return std::sqrt(Constants::epsilon_0 * T_e / n_e / Constants::elementary_charge);
}

template void write_vector_to_binary_file<double>(const std::vector<double>&, const std::string&, size_t, char, bool);
template void write_vector_to_binary_file<int>(const std::vector<int>&, const std::string&, size_t, char, bool);
