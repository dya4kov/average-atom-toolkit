cmake_minimum_required(VERSION 2.8.3)
project (average-atom-tools)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_BINARY_DIR}")

add_subdirectory(numeric-tools)
add_subdirectory(average-atom-tools)
add_subdirectory(tests)
# add_subdirectory (AvAtomPy)
# add_subdirectory (numToolsPy)