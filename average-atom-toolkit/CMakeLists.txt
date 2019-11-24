configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/configuration.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/configuration.h
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/configuration.h
        DESTINATION include/average-atom-toolkit
)
add_subdirectory(thomas-fermi)
add_subdirectory(util)
set(AVERAGE_ATOM_LIBS ${THOMAS_FERMI_LIBS} PARENT_SCOPE)