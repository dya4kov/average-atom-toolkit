add_subdirectory(complete)
set(FERMI_DIRAC_LIBS ${FERMI_DIRAC_COMPLETE_LIB} PARENT_SCOPE)
install(
    FILES complete.h
    DESTINATION include/numeric-tools/specfunc/fermi-dirac
)