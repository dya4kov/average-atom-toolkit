PYTHON_ADD_MODULE(_PyWaveFunction _PyWaveFunction.cxx)
target_link_libraries(_PyWaveFunction
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)

PYTHON_ADD_MODULE(_PyElectronDensityShell _PyElectronDensity.cxx)
target_link_libraries(_PyElectronDensityShell
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)

install(TARGETS _PyWaveFunction
	            _PyElectronDensityShell
        DESTINATION python/pyaatk/TF/atom/shell
)

install (FILES __init__.py DESTINATION python/pyaatk/TF/atom/shell)