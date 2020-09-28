add_library(_PyChemicalPotential MODULE _PyChemicalPotential.cxx)
target_link_libraries(_PyChemicalPotential PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    numtoolslibs
)
set_target_properties(_PyChemicalPotential PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyFreeEnergy MODULE _PyFreeEnergy.cxx)
target_link_libraries(_PyFreeEnergy PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    numtoolslibs
)
set_target_properties(_PyFreeEnergy PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                               SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyChemicalPotential
                _PyFreeEnergy
        DESTINATION python/pyaatk/TF/eos
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/pyaatk/TF/eos
)

add_subdirectory(qx)
add_subdirectory(shell)