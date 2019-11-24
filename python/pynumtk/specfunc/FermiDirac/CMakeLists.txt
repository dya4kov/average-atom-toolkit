add_library(_PyFermiDirac MODULE _PyFermiDirac.cxx)
target_link_libraries(_PyFermiDirac PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyFermiDirac PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                               SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyFermiDirac
        DESTINATION python/pynumtk/specfunc/FermiDirac
)
install (FILES __init__.py DESTINATION python/pynumtk/specfunc/FermiDirac)
