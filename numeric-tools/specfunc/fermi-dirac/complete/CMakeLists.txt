set(FERMI_DIRAC_COMPLETE_LIB FermiDiracComplete PARENT_SCOPE)
add_library(
    FermiDiracComplete 
    SHARED DMHalf.cxx MHalf.cxx Half.cxx ThreeHalf.cxx
)
install(
	TARGETS FermiDiracComplete
    DESTINATION lib
)