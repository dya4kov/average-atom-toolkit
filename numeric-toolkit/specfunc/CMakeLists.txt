add_subdirectory(bessel)
add_subdirectory(fermi-dirac)

set(SPEC_FUNC_LIBS ${BESSEL_LIBS} ${FERMI_DIRAC_LIBS} PARENT_SCOPE)