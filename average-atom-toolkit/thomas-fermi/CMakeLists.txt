add_subdirectory(atom)
add_subdirectory(quantum-exchange)
add_subdirectory(thermodynamics)
set(THOMAS_FERMI_LIBS
  ${THOMAS_FERMI_ATOM_LIBS}
  ${THOMAS_FERMI_QE_LIBS}
  ${THOMAS_FERMI_THERMODYNAMICS_LIBS}
    PARENT_SCOPE
)
install(FILES atom.h DESTINATION include/average-atom-toolkit/thomas-fermi)