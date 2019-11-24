if(ENABLE_MULTITHREADING)
  find_package(Threads)
endif()

add_library(ThomasFermiChemicalPotentialShell SHARED
    chemical-potential.cxx
)
target_link_libraries(
	ThomasFermiChemicalPotentialShell
	ThomasFermiAtom
	ThomasFermiAtomShell
	${CMAKE_THREAD_LIBS_INIT}
	${NUM_TOOLS_LIBS}
)

add_library(ThomasFermiFreeEnergyShell SHARED
    free-energy.cxx
)
target_link_libraries(
	ThomasFermiFreeEnergyShell
	ThomasFermiChemicalPotential
	ThomasFermiChemicalPotentialShell
	ThomasFermiAtom
	ThomasFermiAtomShell
	${CMAKE_THREAD_LIBS_INIT}
	${NUM_TOOLS_LIBS}
)

set(THOMAS_FERMI_EOS_SHELL_LIBS 
	ThomasFermiChemicalPotentialShell
	ThomasFermiFreeEnergyShell
	PARENT_SCOPE
)

install(TARGETS ThomasFermiChemicalPotentialShell
	            ThomasFermiFreeEnergyShell
	    DESTINATION lib)
install(FILES chemical-potential.h
	          free-energy.h
	    DESTINATION include/average-atom-toolkit/thomas-fermi/eos/shell
)