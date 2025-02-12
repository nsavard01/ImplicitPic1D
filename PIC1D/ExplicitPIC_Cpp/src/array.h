#ifndefARRAY_H
#define ARRAY_H

#include <vector>

template<typename T>
class Array {
private:
    std::vector<size_t> dims;
    size_t total_size;
    std::vector<T> data;

public:
    Array(std::vector<size_t> dimensions);
    T& operator()(size_t i, size_t j);
};

#endif  // MULTIARRAY_H