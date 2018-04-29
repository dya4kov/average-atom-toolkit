add_subdirectory(atom)
add_subdirectory(eos)
set(THOMAS_FERMI_LIBS
  ${THOMAS_FERMI_ATOM_LIBS}
  ${THOMAS_FERMI_EOS_LIBS}
    PARENT_SCOPE
)
install(FILES atom.h DESTINATION include/average-atom-toolkit/thomas-fermi)