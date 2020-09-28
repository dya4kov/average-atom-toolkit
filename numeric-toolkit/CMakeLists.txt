add_subdirectory(interpolation)
add_subdirectory(specfunc)
add_subdirectory(ODE)
# set(NUM_TOOLS_LIBS ${SPEC_FUNC_LIBS} ${INTERPOLATION_LIBS} PARENT_SCOPE)

find_package(GSL REQUIRED)
add_library(numtoolslibs INTERFACE)
target_link_libraries(
	numtoolslibs INTERFACE
	${SPEC_FUNC_LIBS} 
	${INTERPOLATION_LIBS}
	GSL::gsl 
    GSL::gslcblas
)
