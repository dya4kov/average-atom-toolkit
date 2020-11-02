if(ENABLE_MULTITHREADING)
add_subdirectory(multithreading)
endif()
add_subdirectory(thomas_fermi)
add_subdirectory(semiclassic)
add_subdirectory(transport)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DESTINATION python
   FILES_MATCHING PATTERN "__init__.py"
)