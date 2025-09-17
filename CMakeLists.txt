# CMake configuration for the project
#
# Dependencies:
# - Eigen: Make sure you have the Eigen library installed.
# - MPI: Ensure MPI is installed for parallel execution.
#
# To install dependencies:
# - Eigen: Install via package manager or manually (http://eigen.tuxfamily.org/).

cmake_minimum_required(VERSION 3.30)
project(YourProjectName)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 
find_package(Eigen3 REQUIRED CONFIG)

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(YourProjectName main.cpp)
target_link_libraries(YourProjectName PRIVATE Eigen3::Eigen MPI::MPI_CXX)
