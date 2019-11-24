add_library(_PyPotentialQX MODULE _PyPotentialQX.cxx)
target_link_libraries(_PyPotentialQX PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${AVERAGE_ATOM_LIBS} 
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyPotentialQX PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyPotentialQX
        DESTINATION python/pyaatk/TF/atom/qx
)
install (FILES __init__.py DESTINATION python/pyaatk/TF/atom/qx)