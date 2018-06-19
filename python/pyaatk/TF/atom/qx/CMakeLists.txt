PYTHON_ADD_MODULE(_PyPotentialQX _PyPotential.cxx)
target_link_libraries(_PyPotentialQX
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)

install(TARGETS _PyPotentialQX
        DESTINATION python/pyaatk/TF/atom/qx
)

install (FILES __init__.py DESTINATION python/pyaatk/TF/atom/qx)