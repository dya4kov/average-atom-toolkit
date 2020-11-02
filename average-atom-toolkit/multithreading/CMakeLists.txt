# use folder name as library name
get_filename_component(LIB_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME ${PROJECT_SHORT_NAME}-${LIB_DIR_NAME})

add_library(${LIB_NAME} INTERFACE)
find_package(Threads REQUIRED)

target_link_libraries(${LIB_NAME} INTERFACE Threads::Threads)
target_include_directories(${LIB_NAME} INTERFACE 
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
)
install(TARGETS ${LIB_NAME} 
	DESTINATION lib
	     EXPORT ${PROJECT_SHORT_NAME}-targets
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	  DESTINATION include/${PROJECT_NAME}
   FILES_MATCHING PATTERN "*.h"
)
