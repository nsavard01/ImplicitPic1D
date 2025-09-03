cmake_minimum_required(VERSION 3.23.0)

set(CMAKE_CXX_COMPILER mpiicpx)
find_package(IntelDPCPP REQUIRED)
project(SSCPIC_1D)

# Set C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Build type and compiler flags -O3 -qopenmp -qmkl -xHost -fast -ipo -Wall
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -qopenmp -qmkl -ffast-math -xHost -ipo -Wall")

# Include directories
include_directories(include)

# Gather all source files
file(GLOB_RECURSE SRC_FILES src/main.cpp
    src/charged_particle_operators/*.cpp
    src/globals/*.cpp
    src/rand_gen/*.cpp
    src/domain/*.cpp
    src/particles/*.cpp
    src/simulation/*.cpp
    src/ES_solvers/*.cpp
    src/solvers/*.cpp
    src/collisions/*.cpp)

# Create the executable
add_executable(sc_pic1d ${SRC_FILES})
