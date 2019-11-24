add_library(
    Interpolation SHARED 
    ${CMAKE_CURRENT_SOURCE_DIR}/base.cxx 
    ${CMAKE_CURRENT_SOURCE_DIR}/linear.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/polynomial.cxx 
    ${CMAKE_CURRENT_SOURCE_DIR}/spline.cxx 
)
install(
	TARGETS Interpolation
    DESTINATION lib
)

set(INTERPOLATION_LIBS Interpolation PARENT_SCOPE)
install(
    FILES base.h linear.h polynomial.h spline.h  
    DESTINATION include/numeric-toolkit/interpolation
)