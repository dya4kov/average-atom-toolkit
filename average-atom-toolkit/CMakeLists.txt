configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/configuration.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/configuration.h
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/configuration.h
        DESTINATION include/${PROJECT_NAME}
)
if(ENABLE_MULTITHREADING)
add_subdirectory(multithreading)
endif()
add_subdirectory(thomas-fermi)
add_subdirectory(semiclassic)
add_subdirectory(util)
add_subdirectory(transport)

install(EXPORT ${PROJECT_SHORT_NAME}-targets 
   DESTINATION share/cmake/${PROJECT_SHORT_NAME}
)
