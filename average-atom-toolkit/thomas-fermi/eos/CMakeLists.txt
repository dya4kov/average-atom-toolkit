# first we add the executable that generates the table
add_executable(table-gen table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  DEPENDS table-gen
  COMMAND table-gen ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)

if(ENABLE_MULTITHREADING)
  find_package(Threads)
endif()

add_library(ThomasFermiChemicalPotential SHARED
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
    chemical-potential.cxx
)
target_link_libraries(
	ThomasFermiChemicalPotential 
	${CMAKE_THREAD_LIBS_INIT}
	${NUM_TOOLS_LIBS}
)

add_library(ThomasFermiFreeEnergy SHARED
	free-energy.cxx
)
target_link_libraries(
	ThomasFermiFreeEnergy
	ThomasFermiChemicalPotential
	${CMAKE_THREAD_LIBS_INIT} 
	${NUM_TOOLS_LIBS}
)

add_subdirectory(qx)
add_subdirectory(shell)

set(THOMAS_FERMI_EOS_LIBS 
	ThomasFermiChemicalPotential
	ThomasFermiFreeEnergy
	${THOMAS_FERMI_EOS_QX_LIBS}
	${THOMAS_FERMI_EOS_SHELL_LIBS}
	PARENT_SCOPE
)

install(TARGETS ThomasFermiChemicalPotential 
	            ThomasFermiFreeEnergy
	    DESTINATION lib)
install(FILES chemical-potential.h
			  free-energy.h 
	    DESTINATION include/average-atom-toolkit/thomas-fermi/eos
)