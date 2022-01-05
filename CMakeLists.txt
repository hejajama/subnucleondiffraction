# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(amplitudelib CXX C)

set(CMAKE_BUILD_TYPE Release)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)
git_local_changes(GIT_LOCAL_CHANGES)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})


# Directory where you have build AmplitudeLib libraries 
link_directories("../amplitudelib_v2/build/lib/")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# AmplitudeLib headers
include_directories("../amplitudelib_v2/")
include_directories("./src/")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Paths to LCPT dipole amplitude code 
# Uncomment and add correct dirs if use, also modify src/ipsat_proton.hpp
#include_directories("/Users/hejajama/code/rhorho")
#link_directories("/Users/hejajama/code/rhorho/build/dipole_interpolation")

#add_subdirectory(libColorDipole/src/)
add_subdirectory(src)
add_subdirectory(tools)

set_target_properties(subnucleondiffraction PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set_target_properties(ddis PROPERTIES
     LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
set_target_properties(mzipsat PROPERTIES
     LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
