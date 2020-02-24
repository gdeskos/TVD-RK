cmake_minimum_required(VERSION 2.8.11)

project(AMR_OFFSHORE_UTILS CXX Fortran)

find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


add_executable(level_set_test
               level_set_test.cpp)
