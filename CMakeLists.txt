cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(AMR_GVOF CXX Fortran)

########################## OPTIONS #####################################
#Options for the executable
option(AMR_WIND_ENABLE_MPI "Enable MPI" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

########################### AMReX #####################################

set(AMREX_SUBMOD_LOCATION "submods/amrex_new")
add_subdirectory(${AMREX_SUBMOD_LOCATION})

if(AMR_WIND_ENABLE_FCOMPARE OR AMR_WIND_TEST_WITH_FCOMPARE)
      add_subdirectory(${AMREX_SUBMOD_LOCATION}/Tools/Plotfile)
endif()

find_package(MPI REQUIRED)

#Options for C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# General information about machine, compiler, and build type
message(STATUS "amr-gvof Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(AMREX_SUBMOD_LOCATION "submods/amrex_new")

set(CMAKE_CXX_COMPILER ${AMREX_CXX_COMPILER})

set(amr_gvof_exe_name "amr_gvof")
add_executable(${amr_gvof_exe_name}
               src/main.cpp
               src/AmrCoreAdv.cpp
               src/AdvancePhiAllLevels.cpp
               src/AdvancePhiAtLevel.cpp
               src/DefineVelocity.cpp)
target_link_libraries(${amr_gvof_exe_name} amrex)

