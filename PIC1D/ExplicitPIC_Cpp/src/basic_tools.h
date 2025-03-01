#ifndef BASIC_TOOLS_H
#define BASIC_TOOLS_H

#include <string>
#include <vector>

template <typename T>
void write_vector_to_binary_file(const std::vector<T>& vec, const std::string& filename, size_t paddingBytes = 0, char paddingValue = 0x00);

bool createDirectory(const std::string& dirName);

bool directoryExists(const std::string& dirName);

#endif 