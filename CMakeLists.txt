# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(amplitudelib CXX C Fortran)

set(CMAKE_BUILD_TYPE Release)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)
git_local_changes(GIT_LOCAL_CHANGES)

# Directory where you have build AmplitudeLib libraries 
link_directories("../amplitudelib_v2/build/lib/")

# AmplitudeLib headers
include_directories("../amplitudelib_v2/")

add_subdirectory(libColorDipole/src/)
add_subdirectory(src)


set_target_properties(subnucleondiffraction PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )