add_subdirectory(atom)
set(THOMAS_FERMI_LIBS ${THOMAS_FERMI_ATOM_LIBS} PARENT_SCOPE)
install (FILES atom.h DESTINATION include/average-atom-tools/thomas-fermi)