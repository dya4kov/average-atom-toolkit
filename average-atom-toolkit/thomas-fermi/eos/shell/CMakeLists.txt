# add_subdirectory(ODE)
set(EOS_SHELL_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical-potential.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/free-energy.cxx
    PARENT_SCOPE
)
# install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h
#   DESTINATION include/${PROJECT_NAME}/${LIB_DIR_NAME}/eos/shell
# )