add_library(
    FermiDirac SHARED 
    ${CMAKE_CURRENT_SOURCE_DIR}/complete/DMHalf.cxx 
    ${CMAKE_CURRENT_SOURCE_DIR}/complete/MHalf.cxx 
    ${CMAKE_CURRENT_SOURCE_DIR}/complete/Half.cxx 
    ${CMAKE_CURRENT_SOURCE_DIR}/complete/ThreeHalf.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/Yfunction/Yfunction.cxx
)
install(
	TARGETS FermiDirac
    DESTINATION lib
)

set(FERMI_DIRAC_LIBS FermiDirac PARENT_SCOPE)
install(
    FILES complete.h Yfunction.h
    DESTINATION include/numeric-tools/specfunc/fermi-dirac
)