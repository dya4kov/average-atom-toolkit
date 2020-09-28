add_library(_PyChemicalPotentialShell MODULE _PyChemicalPotentialShell.cxx)
target_link_libraries(_PyChemicalPotentialShell PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    numtoolslibs
)
set_target_properties(_PyChemicalPotentialShell PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                           SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyFreeEnergyShell MODULE _PyFreeEnergyShell.cxx)
target_link_libraries(_PyFreeEnergyShell PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    numtoolslibs
)
set_target_properties(_PyFreeEnergyShell PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                    SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyChemicalPotentialShell
                _PyFreeEnergyShell
        DESTINATION python/pyaatk/TF/eos/shell
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/pyaatk/TF/eos/shell
)