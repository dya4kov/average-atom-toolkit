# find Threads library
find_package(Threads)
# add ThomasFermiAtom atom library
add_library(ThomasFermiThermodynamics SHARED
	free-energy.cxx
	qe-free-energy.cxx
)
target_link_libraries(
	ThomasFermiThermodynamics
	ThomasFermiAtom
	ThomasFermiQuantumExchange
	${CMAKE_THREAD_LIBS_INIT} 
	${NUM_TOOLS_LIBS}
)

set(THOMAS_FERMI_THERMODYNAMICS_LIBS ThomasFermiThermodynamics PARENT_SCOPE)

install(TARGETS ThomasFermiThermodynamics DESTINATION lib)
install(FILES free-energy.h qe-free-energy.h
	    DESTINATION include/average-atom-tools/thomas-fermi/thermodynamics)