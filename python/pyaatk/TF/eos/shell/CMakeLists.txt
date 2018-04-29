PYTHON_ADD_MODULE(_PyChemicalPotentialShell _PyChemicalPotentialShell.cxx)
target_link_libraries(_PyChemicalPotentialShell
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyFreeEnergyShell _PyFreeEnergyShell.cxx)
target_link_libraries(_PyFreeEnergyShell
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
install(TARGETS _PyChemicalPotentialShell
                _PyFreeEnergyShell
        DESTINATION python/pyaatk/TF/eos/shell
)

install (FILES __init__.py DESTINATION python/pyaatk/TF/eos/shell)