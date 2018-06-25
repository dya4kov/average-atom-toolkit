if(ENABLE_MULTITHREADING)
  find_package(Threads)
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical-potential.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/chemical-potential.h
)
add_library(ThomasFermiChemicalPotentialShell SHARED
    chemical-potential.cxx
)
target_link_libraries(
	ThomasFermiChemicalPotentialShell
	ThomasFermiAtom
	${CMAKE_THREAD_LIBS_INIT}
	${NUM_TOOLS_LIBS}
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/free-energy.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/free-energy.h
)
add_library(ThomasFermiFreeEnergyShell SHARED
    free-energy.cxx
)
target_link_libraries(
	ThomasFermiFreeEnergyShell
	ThomasFermiChemicalPotential
	ThomasFermiChemicalPotentialShell
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
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/chemical-potential.h
	          ${CMAKE_CURRENT_BINARY_DIR}/free-energy.h
	    DESTINATION include/average-atom-toolkit/thomas-fermi/eos/shell
)