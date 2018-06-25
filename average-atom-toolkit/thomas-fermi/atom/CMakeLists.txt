if(ENABLE_MULTITHREADING)
  # find Threads library
  find_package(Threads)
endif()
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/energy-level.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/energy-level.h
)
# add ThomasFermiAtom atom library
add_library(ThomasFermiAtom SHARED
    potential.cxx
    rotate-points.cxx
    action.cxx
    energy-level.cxx
    electron-density.cxx
    electron-states.cxx
    atom.cxx
)
target_link_libraries(
  ThomasFermiAtom
  ThomasFermiChemicalPotential
  ${CMAKE_THREAD_LIBS_INIT}
  ${NUM_TOOLS_LIBS}
)
add_subdirectory(qx)
add_subdirectory(shell)
set(THOMAS_FERMI_ATOM_LIBS ThomasFermiAtom ThomasFermiAtomQX ThomasFermiAtomShell PARENT_SCOPE)

install(TARGETS ThomasFermiAtom DESTINATION lib)
install(FILES action.h
              electron-density.h 
              electron-states.h
              ${CMAKE_CURRENT_BINARY_DIR}/energy-level.h 
              potential.h 
              rotate-points.h
        DESTINATION include/average-atom-toolkit/thomas-fermi/atom
)