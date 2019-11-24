# first we add the executable that generates the table
add_executable(table-gen-qx table-gen.cxx)
# add the command to generate the source code
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
  DEPENDS table-gen-qx
  COMMAND table-gen-qx ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)

if(ENABLE_MULTITHREADING)
  find_package(Threads)
endif()

add_library(ThomasFermiChemicalPotentialQX SHARED
  ${CMAKE_CURRENT_BINARY_DIR}/table.cxx
    chemical-potential.cxx
)
target_link_libraries(
	ThomasFermiChemicalPotentialQX 
	ThomasFermiChemicalPotential
	${CMAKE_THREAD_LIBS_INIT}
	${NUM_TOOLS_LIBS}
)

add_library(ThomasFermiFreeEnergyQX SHARED
    free-energy.cxx
)
target_link_libraries(
	ThomasFermiFreeEnergyQX
	ThomasFermiChemicalPotential
	ThomasFermiChemicalPotentialQX
	${CMAKE_THREAD_LIBS_INIT}
	${NUM_TOOLS_LIBS}
)

set(THOMAS_FERMI_EOS_QX_LIBS 
	ThomasFermiChemicalPotentialQX
	ThomasFermiFreeEnergyQX
	PARENT_SCOPE
)

install(TARGETS ThomasFermiChemicalPotentialQX
	            ThomasFermiFreeEnergyQX
	    DESTINATION lib)
install(FILES chemical-potential.h
	          free-energy.h
	    DESTINATION include/average-atom-toolkit/thomas-fermi/eos/qx
)