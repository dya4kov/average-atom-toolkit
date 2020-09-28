# add semiclassic atom library
add_subdirectory(potential)
add_library(aatk-transport SHARED
    ${AATK_TRANSPORT_POTENTIAL_SOURCES}
    delta.cxx
    cross-section.cxx
    Kfunction.cxx
)
target_link_libraries(
  aatk-transport PUBLIC
  numtoolslibs
)

install(TARGETS aatk-transport DESTINATION lib)
