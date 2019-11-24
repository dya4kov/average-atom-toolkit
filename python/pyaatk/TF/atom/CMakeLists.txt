add_library(_PyPotential MODULE _PyPotential.cxx)
target_link_libraries(_PyPotential PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyPotential PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                              SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyElectronDensity MODULE _PyElectronDensity.cxx)
target_link_libraries(_PyElectronDensity PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyElectronDensity PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                    SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyElectronStates MODULE _PyElectronStates.cxx)
target_link_libraries(_PyElectronStates PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyElectronStates PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                   SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyEnergyLevel MODULE _PyEnergyLevel.cxx)
target_link_libraries(_PyEnergyLevel PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyEnergyLevel PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyRotatePoints MODULE _PyRotatePoints.cxx)
target_link_libraries(_PyRotatePoints PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyRotatePoints PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                 SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyPotential
                _PyElectronDensity 
                _PyRotatePoints
                _PyEnergyLevel
                _PyElectronStates
        DESTINATION python/pyaatk/TF/atom
)
install (FILES __init__.py DESTINATION python/pyaatk/TF/atom)
add_subdirectory(qx)
add_subdirectory(shell)