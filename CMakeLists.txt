cmake_minimum_required(VERSION 3.15.0)
project (average-atom-toolkit)

set(CMAKE_CXX_STANDARD             11)
set(PROJECT_SHORT_NAME           aatk)
set(CMAKE_CXX_FLAGS_RELEASE     "-O3")

option(CMAKE_BUILD_TYPE         Debug)
option(ENABLE_MULTITHREADING       ON)
option(ENABLE_PYTHON               ON)
option(ENABLE_TESTS                ON)

add_subdirectory(average-atom-toolkit)

if(ENABLE_TESTS)
add_subdirectory(tests)
endif()

if(ENABLE_PYTHON)
  set(PY_PROJECT_NAME pyaatk)
  add_subdirectory(python)
endif()

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/aatkConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/aatkConfig.cmake"
  INSTALL_DESTINATION "share/cmake/aatk"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install the configuration file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/aatkConfig.cmake
  DESTINATION share/cmake/aatk
)