add_library(gsl_sf_bessel SHARED IMPORTED GLOBAL)
if (NOT APPLE)
set_target_properties(gsl_sf_bessel PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libgsl_sf_bessel.so)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/lib/libgsl_sf_bessel.so DESTINATION lib)
else()
set_target_properties(gsl_sf_bessel PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libgsl_sf_bessel.dylib)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/lib/libgsl_sf_bessel.dylib DESTINATION lib)
endif()
set(BESSEL_LIBS gsl_sf_bessel PARENT_SCOPE)
install(FILES I0.h I1.h In.h Inu.h
              J0.h J1.h Jn.h Jnu.h
              K0.h K1.h Kn.h Knu.h
              Y0.h Y1.h Yn.h Ynu.h
  DESTINATION include/numeric-toolkit/specfunc/bessel
)