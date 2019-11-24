add_subdirectory(interpolation)
add_subdirectory(specfunc)
add_subdirectory(ODE)
set(NUM_TOOLS_LIBS ${SPEC_FUNC_LIBS} ${INTERPOLATION_LIBS} PARENT_SCOPE)
