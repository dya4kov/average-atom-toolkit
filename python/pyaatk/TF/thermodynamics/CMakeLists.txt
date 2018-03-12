PYTHON_ADD_MODULE(_PyChemicalPotential _PyChemicalPotential.cxx)
target_link_libraries(_PyChemicalPotential
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyFreeEnergy _PyFreeEnergy.cxx)
target_link_libraries(_PyFreeEnergy
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
install(TARGETS _PyChemicalPotential
                _PyFreeEnergy
        DESTINATION python/pyaatk/TF/thermodynamics
)
install (FILES __init__.py DESTINATION python/pyaatk/TF/thermodynamics)
add_subdirectory(QE)
add_subdirectory(shell)