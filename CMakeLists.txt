# Copyright CEA (2024)
# Contributors: Théophile Bonnet, Davide Mancusi
# tlfab2@cam.ac.uk, davide.mancusi@cea.fr
# 
# This software is a computer program whose purpose is to investigate
# correlations in neutron transport calculations in very simplified
# conditions.
# 
# This software is governed by the CeCILL  license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

cmake_minimum_required(VERSION 3.11)
project(coyote CXX)

option(USE_OMP "Compile with OpenMP for shared memory parralelisation" ON)
include(FetchContent)

set(COYOTE_SOURCES 
    sources/simulation_class.cpp 
    sources/ressources.cpp 
    sources/tally.cpp 
    sources/particle_class.cpp 
    sources/scores_class.cpp
    sources/geometry_class.cpp
    sources/region_class.cpp
    sources/population_control_class.cpp
    sources/collision_kernel_class.cpp
    sources/decay_kernel_class.cpp
    sources/vector_class.cpp
    sources/front_screen_class.cpp
    sources/parser.cpp
    sources/settings.cpp
    sources/material_class.cpp
    sources/random_class.cpp
    sources/source_class.cpp
    sources/error.cpp
    sources/static/static_collision_kernel_class.cpp
    sources/static/static_scores_class.cpp
    sources/static/static_simulation_class.cpp
    sources/static/static_population_control.cpp)

# Get YAML-CPP version 0.7.0
message(STATUS "Downloading yaml-cpp 0.7.0")
FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG        yaml-cpp-0.7.0)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Enable testing")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Enable parse tools")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Enable contrib stuff in library")
set(YAML_CPP_INSTALL OFF CACHE BOOL "Enable generation of install target")
FetchContent_MakeAvailable(yaml-cpp)

# Get NDArray
message(STATUS "Downloading NDArray")
FetchContent_Declare(NDArray
  GIT_REPOSITORY https://github.com/HunterBelanger/ndarray.git
  GIT_TAG        master
)
set(NDARRAY_INSTALL OFF CACHE BOOL "Install NDArray")
FetchContent_MakeAvailable(NDArray)

# Get PCG-C++
message(STATUS "Downloading PCG-C++")
FetchContent_Declare(PCG_CXX
  GIT_REPOSITORY https://github.com/HunterBelanger/pcg-cpp.git
  GIT_TAG        feature/cmake
)
set(PCG_CPP_INSTALL OFF CACHE BOOL "Install PCG-CPP")
FetchContent_MakeAvailable(PCG_CXX)

# Get the Git hash for define
add_executable(coyote sources/main.cpp ${COYOTE_SOURCES})
target_include_directories(coyote PRIVATE include)
target_compile_features(coyote PRIVATE cxx_std_17)
target_compile_options(coyote PRIVATE -W -Wall -Wextra -Wpedantic -Wno-unused-but-set-parameter)
target_link_libraries(coyote PUBLIC yaml-cpp NDArray::NDArray PCG::PCG_CXX)

# Find OpenMP if desired
#if(USE_OMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(coyote PUBLIC OpenMP::OpenMP_CXX)
    endif()
#endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Allow for install of executable
install(TARGETS coyote
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

