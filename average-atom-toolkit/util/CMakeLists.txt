# use folder name as library name
get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME "${PROJECT_SHORT_NAME}-${DIR_NAME}")
add_subdirectory(specfunc)
add_library(${LIB_NAME} SHARED ${SPECFUNC_SRC})
target_link_libraries(${LIB_NAME} PUBLIC
	GSL::gsl 
    GSL::gslcblas
)
target_include_directories(${LIB_NAME} PUBLIC 
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)
install(TARGETS ${LIB_NAME}
	DESTINATION lib
         EXPORT ${PROJECT_SHORT_NAME}-targets
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DESTINATION include/${PROJECT_NAME}
   FILES_MATCHING PATTERN "*.h"
)