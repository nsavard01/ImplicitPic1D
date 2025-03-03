#include "basic_tools.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>

template <typename T>
void write_vector_to_binary_file(const std::vector<T>& vec, const std::string& filename, size_t paddingBytes, char paddingValue) {
    // Calculate total size (data + padding)
    size_t size_of_T = sizeof(T);
    size_t total_size = paddingBytes + vec.size() * size_of_T;
    

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

    // Step 5: Write the byte buffer to the binary file
    std::ofstream outFile(filename, std::ios::binary);
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

template void write_vector_to_binary_file<double>(const std::vector<double>&, const std::string&, size_t, char);
template void write_vector_to_binary_file<int>(const std::vector<int>&, const std::string&, size_t, char);
