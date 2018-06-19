PYTHON_ADD_MODULE(_PyPotential _PyPotential.cxx)
target_link_libraries(_PyPotential 
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyElectronDensity _PyElectronDensity.cxx)
target_link_libraries(_PyElectronDensity
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyRotatePoints _PyRotatePoints.cxx)
target_link_libraries(_PyRotatePoints
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyEnergyLevel _PyEnergyLevel.cxx)
target_link_libraries(_PyEnergyLevel
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyElectronStates _PyElectronStates.cxx)
target_link_libraries(_PyElectronStates
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
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