# use folder name as library name
get_filename_component(LIB_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME ${PROJECT_SHORT_NAME}-${LIB_DIR_NAME})
# gather sources
add_subdirectory(thomas-fermi)
add_subdirectory(semiclassic)
# add library
add_library(${LIB_NAME} 
	SHARED base.cxx 
	     ${ATOM_TF_SRC} 
	     ${ATOM_SC_SRC}
)
add_dependencies(${LIB_NAME} tftable-gen tftablegenerated)
target_link_libraries(${LIB_NAME} aatk-util)
target_include_directories(${LIB_NAME} PUBLIC 
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
)
install(TARGETS ${LIB_NAME} 
	DESTINATION lib
	     EXPORT ${PROJECT_SHORT_NAME}-targets
)
# copy headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DESTINATION include/${PROJECT_NAME}
   FILES_MATCHING PATTERN "*.h"
)