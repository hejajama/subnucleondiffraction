# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# project name and enable C++ support
project(subnucleondiffraction CXX C)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)

# Required to compile on MacOS with homebrew installed OpenMP
cmake_policy(SET CMP0074 NEW)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)
git_local_changes(GIT_LOCAL_CHANGES)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})


find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Paths to LCPT dipole amplitude code 
# Uncomment and add correct dirs if use, also modify src/ipsat_proton.hpp and src/CMakeLists.txt
#include_directories("/Users/hejajama/code/rhorho")
#link_directories("/Users/hejajama/code/rhorho/build/dipole_interpolation")

# IPsat fit from arXiv:1212.2974
# not compiled by default
#add_subdirectory(libColorDipole/src/)

add_subdirectory(src)
add_subdirectory(tools)

set_target_properties(subnucleondiffraction PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set_target_properties(test PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

set_target_properties(ddis PROPERTIES
     LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
set_target_properties(mzipsat PROPERTIES
     LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
