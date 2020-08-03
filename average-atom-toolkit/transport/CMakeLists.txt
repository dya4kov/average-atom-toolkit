# add semiclassic atom library
find_package(GSL REQUIRED)    # See below (2)
add_subdirectory(potential)
add_library(aatk-transport SHARED
    ${AATK_TRANSPORT_POTENTIAL_SOURCES}
    delta.cxx
    cross-section.cxx
    Kfunction.cxx
)
target_link_libraries(
  aatk-transport PUBLIC
  ${NUM_TOOLS_LIBS}
  GSL::gsl 
  GSL::gslcblas
)

install(TARGETS aatk-transport DESTINATION lib)
