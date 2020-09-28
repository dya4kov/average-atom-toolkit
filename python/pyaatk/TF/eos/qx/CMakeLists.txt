add_library(_PyChemicalPotentialQX MODULE _PyChemicalPotentialQX.cxx)
target_link_libraries(_PyChemicalPotentialQX PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    numtoolslibs
)
set_target_properties(_PyChemicalPotentialQX PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyFreeEnergyQX MODULE _PyFreeEnergyQX.cxx)
target_link_libraries(_PyFreeEnergyQX PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    numtoolslibs
)
set_target_properties(_PyFreeEnergyQX PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                 SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyChemicalPotentialQX
                _PyFreeEnergyQX
        DESTINATION python/pyaatk/TF/eos/qx
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/pyaatk/TF/eos/qx
)