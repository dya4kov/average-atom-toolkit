# add_subdirectory(ODE)
set(EOS_SHELL_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical-potential.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/free-energy.cxx
    PARENT_SCOPE
)