# dummy config to create file in build directory
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/table.cxx)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/table.cxx.in
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
)
endif()
# first we add the executable that generates the table
add_executable(table-gen-qx table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  COMMAND table-gen-qx ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS table-gen-qx ${CMAKE_CURRENT_SOURCE_DIR}/table.cxx.in
)
add_custom_target(
    tableqxgenerated ALL
    DEPENDS table-gen-qx ${CMAKE_CURRENT_BINARY_DIR}/table.cxx 
)
# add_subdirectory(ODE)
set(EOS_QX_SRC
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical-potential.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/free-energy.cxx
    PARENT_SCOPE
)
