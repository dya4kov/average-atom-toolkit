PYTHON_ADD_MODULE(_PyChemicalPotentialQX _PyChemicalPotentialQX.cxx)
target_link_libraries(_PyChemicalPotentialQX
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyFreeEnergyQX _PyFreeEnergyQX.cxx)
target_link_libraries(_PyFreeEnergyQX
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
install(TARGETS _PyChemicalPotentialQX
                _PyFreeEnergyQX
        DESTINATION python/pyaatk/TF/eos/qx
)

install (FILES __init__.py DESTINATION python/pyaatk/TF/eos/qx)