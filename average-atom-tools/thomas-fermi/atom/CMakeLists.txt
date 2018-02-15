# first we add the executable that generates the table
add_executable(table-gen table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  DEPENDS table-gen
  COMMAND table-gen ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)
# add ODE headers
add_subdirectory(ODE)
# find Threads library
find_package(Threads)
# add ThomasFermiAtom atom library
add_library(ThomasFermiAtom SHARED
    ${CMAKE_CURRENT_BINARY_DIR}/table.cxx 
    potential.cxx
    rotate-points.cxx
    action.cxx
    energy-level.cxx
    electron-density.cxx
    electron-states.cxx
    atom.cxx
)
target_link_libraries(ThomasFermiAtom ${CMAKE_THREAD_LIBS_INIT} ${NUM_TOOLS_LIBS})

set(THOMAS_FERMI_ATOM_LIBS ThomasFermiAtom PARENT_SCOPE)

install(TARGETS ThomasFermiAtom DESTINATION lib)
install(FILES action.h
              electron-states.h 
              energy-level.h 
              potential.h 
              rotate-points.h
              table.h 
        DESTINATION include/average-atom-tools/thomas-fermi/atom
)