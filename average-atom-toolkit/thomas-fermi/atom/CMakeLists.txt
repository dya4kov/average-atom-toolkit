# find Threads library
find_package(Threads)
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

set(THOMAS_FERMI_ATOM_LIBS ThomasFermiAtom PARENT_SCOPE)

install(TARGETS ThomasFermiAtom DESTINATION lib)
install(FILES action.h
              electron-density.h 
              electron-states.h
              energy-level.h 
              potential.h 
              rotate-points.h
        DESTINATION include/average-atom-toolkit/thomas-fermi/atom
)