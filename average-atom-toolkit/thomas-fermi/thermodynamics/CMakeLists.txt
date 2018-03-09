# first we add the executable that generates the table
add_executable(table-gen table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  DEPENDS table-gen
  COMMAND table-gen ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)
# find Threads library
find_package(Threads)
add_library(ThomasFermiChemicalPotential SHARED
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
    chemical-potential.cxx
)
target_link_libraries(
	ThomasFermiChemicalPotential 
	${CMAKE_THREAD_LIBS_INIT}
	${NUM_TOOLS_LIBS}
)
# add ThomasFermiAtom atom library
add_library(ThomasFermiThermodynamics SHARED
	free-energy.cxx
	qe-free-energy.cxx
)
target_link_libraries(
	ThomasFermiThermodynamics
	ThomasFermiChemicalPotential
	ThomasFermiQuantumExchange
	${CMAKE_THREAD_LIBS_INIT} 
	${NUM_TOOLS_LIBS}
)

set(THOMAS_FERMI_THERMODYNAMICS_LIBS 
	ThomasFermiChemicalPotential 
	ThomasFermiThermodynamics 
	PARENT_SCOPE
)

install(TARGETS ThomasFermiThermodynamics 
	            ThomasFermiChemicalPotential 
	    DESTINATION lib)
install(FILES free-energy.h 
	          qe-free-energy.h
	    DESTINATION include/average-atom-toolkit/thomas-fermi/thermodynamics
)