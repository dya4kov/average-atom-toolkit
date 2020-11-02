set(ATOM_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/potential.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/rotate-points.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/action.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/energy-level.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/electron-density.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/electron-states.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/atom.cxx
)
add_subdirectory(qx)
add_subdirectory(shell)
# gather all sources and pass them to parent CMakeLists.txt
set(ATOM_SRC
  ${ATOM_SRC}
  ${ATOM_QX_SRC}
  ${ATOM_SHELL_SRC}
    PARENT_SCOPE
)