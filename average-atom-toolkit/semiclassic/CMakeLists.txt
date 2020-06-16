# add semiclassic atom library
add_library(aatk-semiclassic-atom SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/potential.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/rotate-points.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/action.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/energy-levels.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/electron-density.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/electron-states.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/wave-function.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/atom/atom.cxx
)
target_link_libraries(
  aatk-semiclassic-atom PUBLIC
  ThomasFermiAtom
  ThomasFermiChemicalPotential
  ${NUM_TOOLS_LIBS}
)
if(ENABLE_MULTITHREADING)
  target_link_libraries(aatk-semiclassic-atom PUBLIC aatk-multithreading)
endif()

install(TARGETS aatk-semiclassic-atom DESTINATION lib)
install(FILES atom.h
        DESTINATION include/average-atom-toolkit/semiclassic
)