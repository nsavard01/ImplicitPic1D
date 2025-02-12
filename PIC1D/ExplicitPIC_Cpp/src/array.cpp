
#include "array.h"

template<typename T>
Array<T>::Array(std::vector<size_t> dimensions) : dims(dimensions) {
    total_size = 1;
    for (size_t d : dims) total_size *= d;
    data.resize(total_size);
}

template<typename T>
T& Array<T>::operator()(size_t i, size_t j) {
    return data[i * dims[1] + j];  // Row-major access
}
