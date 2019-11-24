add_library(_PyPotentialShell MODULE _PyPotential.cxx)
target_link_libraries(_PyPotentialShell PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyPotentialShell PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                   SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyWaveFunction MODULE _PyWaveFunction.cxx)
target_link_libraries(_PyWaveFunction PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyWaveFunction PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                 SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(_PyElectronDensityShell MODULE _PyElectronDensity.cxx)
target_link_libraries(_PyElectronDensityShell PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyElectronDensityShell PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyWaveFunction
	            _PyElectronDensityShell
                _PyPotentialShell
        DESTINATION python/pyaatk/TF/atom/shell
)

install (FILES __init__.py DESTINATION python/pyaatk/TF/atom/shell)