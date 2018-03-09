# first we add the executable that generates the table
add_executable(qe-table-gen table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  DEPENDS qe-table-gen
  COMMAND qe-table-gen ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)
# add ODE headers
add_subdirectory(ODE)
# find Threads library
# find_package(Threads)
# add ThomasFermiQuantumExchange atom library
add_library(ThomasFermiQuantumExchange SHARED
	potential.cxx
	${CMAKE_CURRENT_BINARY_DIR}/table.cxx
)
target_link_libraries(
	ThomasFermiQuantumExchange 
	ThomasFermiChemicalPotential 
	${CMAKE_THREAD_LIBS_INIT} 
	${NUM_TOOLS_LIBS}
)

set(THOMAS_FERMI_QE_LIBS ThomasFermiQuantumExchange PARENT_SCOPE)

install(TARGETS ThomasFermiQuantumExchange DESTINATION lib)
install(FILES table.h potential.h
	    DESTINATION include/average-atom-toolkit/thomas-fermi/quantum-exchange)