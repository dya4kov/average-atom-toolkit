# use folder name as library name
get_filename_component(LIB_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME ${PROJECT_SHORT_NAME}-${LIB_DIR_NAME})
# find numeric-toolkit
find_package(numtk)
# gather sources
add_subdirectory(atom)
add_subdirectory(eos)
# add library
add_library(${LIB_NAME} SHARED ${ATOM_SRC} ${EOS_SRC})
add_dependencies(${LIB_NAME} table-gen table-gen-qx tablegenerated tableqxgenerated )
target_link_libraries(${LIB_NAME} 
	PUBLIC numtk-specfunc
	       numtk-interpolation 
	INTERFACE numtk-ODE
)
if(ENABLE_MULTITHREADING)
  # find Threads library
  find_package(Threads)
  # link to pthreads
  target_link_libraries(${LIB_NAME} PUBLIC Threads::Threads)
endif()
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
