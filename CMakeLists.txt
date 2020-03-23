# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(amplitudelib CXX C) # Fortran)

set(CMAKE_BUILD_TYPE Release)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)
git_local_changes(GIT_LOCAL_CHANGES)

find_package(GSL REQUIRED) 

# Directory where you have build AmplitudeLib libraries 
link_directories("../amplitudelib_oberon/build/lib/")

# AmplitudeLib headers
include_directories("../amplitudelib_v2/")
include_directories("./src/")

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
