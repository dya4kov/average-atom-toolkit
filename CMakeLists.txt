cmake_minimum_required(VERSION 3.1.0)
project (average-atom-toolkit)

set(CMAKE_CXX_STANDARD            11)
set(CMAKE_BUILD_TYPE         Release)
set(CMAKE_CXX_FLAGS_RELEASE    "-O3")
set(ENABLE_MULTITHREADING        OFF)
set(ENABLE_PYTHON                 ON)

include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_BINARY_DIR}")

add_subdirectory(numeric-toolkit)
add_subdirectory(average-atom-toolkit)
add_subdirectory(tests)

if(ENABLE_PYTHON)
  add_subdirectory(python)
endif()