PYTHON_ADD_MODULE(_PyFD3half _PyFD3half.cxx)
target_link_libraries(_PyFD3half 
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyFDhalf _PyFDhalf.cxx)
target_link_libraries(_PyFDhalf
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyFDmhalf _PyFDmhalf.cxx)
target_link_libraries(_PyFDmhalf
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${NUM_TOOLS_LIBS}
)
PYTHON_ADD_MODULE(_PyFDdmhalf _PyFDdmhalf.cxx)
target_link_libraries(_PyFDdmhalf
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} 
    ${NUM_TOOLS_LIBS}
)

install(TARGETS _PyFD3half
                _PyFDhalf 
                _PyFDmhalf
                _PyFDdmhalf
        DESTINATION python/pynumtk/specfunc/FermiDirac
)
install (FILES __init__.py DESTINATION python/pynumtk/specfunc/FermiDirac)
