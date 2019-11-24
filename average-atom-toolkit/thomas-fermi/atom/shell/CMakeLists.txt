# find Threads library
find_package(Threads)
#
add_library(ThomasFermiAtomShell SHARED
    wave-function.cxx
    electron-density.cxx
)
target_link_libraries(
  ThomasFermiAtomShell
  ThomasFermiAtom
  ThomasFermiChemicalPotential
  ${CMAKE_THREAD_LIBS_INIT}
  ${NUM_TOOLS_LIBS}
)
add_library(ThomasFermiAtomShellPotential SHARED
    potential.cxx
)
target_link_libraries(
  ThomasFermiAtomShellPotential
  ThomasFermiAtom
  ThomasFermiAtomShell
  ThomasFermiChemicalPotential
  ThomasFermiChemicalPotentialShell
  ${CMAKE_THREAD_LIBS_INIT}
  ${NUM_TOOLS_LIBS}
)
install(TARGETS ThomasFermiAtomShell ThomasFermiAtomShellPotential DESTINATION lib)
install(FILES wave-function.h
	            electron-density.h
	            potential.h
        DESTINATION include/average-atom-toolkit/thomas-fermi/atom/shell
)