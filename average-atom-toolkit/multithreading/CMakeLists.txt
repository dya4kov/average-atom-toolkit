add_library(aatk-multithreading INTERFACE)
find_package(Threads REQUIRED)

target_link_libraries(aatk-multithreading INTERFACE Threads::Threads)
target_include_directories(aatk-multithreading
    INTERFACE 
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
)

install(FILES 
        ${CMAKE_CURRENT_SOURCE_DIR}/thread-pool.h
	    DESTINATION include/average-atom-toolkit/multithreading
)
