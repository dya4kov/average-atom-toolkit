# dummy config to create file in build directory
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/table.cxx)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/table.cxx.in
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
)
endif()
# first we add the executable that generates the table
add_executable(table-gen table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  COMMAND table-gen ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS table-gen ${CMAKE_CURRENT_SOURCE_DIR}/table.cxx.in
)
add_custom_target(
    tablegenerated ALL
    DEPENDS table-gen ${CMAKE_CURRENT_BINARY_DIR}/table.cxx 
)
set(EOS_SRC
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical-potential.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/free-energy.cxx
)

add_subdirectory(qx)
add_subdirectory(shell)
# gather all sources and pass them to parent CMakeLists.txt
set(EOS_SRC 
  ${EOS_SRC}
  ${EOS_QX_SRC}
  ${EOS_SHELL_SRC}
    PARENT_SCOPE
)