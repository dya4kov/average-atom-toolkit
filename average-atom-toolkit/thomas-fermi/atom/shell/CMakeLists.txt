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
install(TARGETS ThomasFermiAtomShell DESTINATION lib)
install(FILES wave-function.h
	          electron-density.h
        DESTINATION include/average-atom-toolkit/thomas-fermi/atom/shell
)