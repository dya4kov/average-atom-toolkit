# dummy config to create file in build directory
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/table.cxx)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/table.cxx.in
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
)
endif()
# first we add the executable that generates the table
add_executable(tftable-gen table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  COMMAND tftable-gen ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS tftable-gen ${CMAKE_CURRENT_SOURCE_DIR}/table.cxx.in
)
add_custom_target(
    tftablegenerated ALL
    DEPENDS tftable-gen ${CMAKE_CURRENT_BINARY_DIR}/table.cxx 
)

set(ATOM_TF_SRC
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/potential.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/interface.cxx
    PARENT_SCOPE
)