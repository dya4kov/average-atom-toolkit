add_library(_PyInterpolation MODULE _PyInterpolation.cxx)
target_link_libraries(_PyInterpolation PRIVATE pybind11::module 
    ${PYTHON_LIBRARIES}
    ${NUM_TOOLS_LIBS}
)
set_target_properties(_PyInterpolation PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                                  SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS _PyInterpolation
        DESTINATION python/pynumtk/interpolation
)
install (FILES __init__.py DESTINATION python/pynumtk/interpolation)
