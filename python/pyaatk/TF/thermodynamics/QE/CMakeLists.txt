PYTHON_ADD_MODULE(_PyChemicalPotentialQE _PyChemicalPotentialQE.cxx)
target_link_libraries(_PyChemicalPotentialQE
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyFreeEnergyQE _PyFreeEnergyQE.cxx)
target_link_libraries(_PyFreeEnergyQE
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
install(TARGETS _PyChemicalPotentialQE
                _PyFreeEnergyQE
        DESTINATION python/pyaatk/TF/thermodynamics/QE
)

install (FILES __init__.py DESTINATION python/pyaatk/TF/thermodynamics/QE)