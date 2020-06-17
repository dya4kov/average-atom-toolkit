add_library(_aatk_multithreading MODULE _multithreading.cxx)
target_link_libraries(_aatk_multithreading 
    PRIVATE pybind11::module 
    PRIVATE aatk-multithreading
)
set_target_properties(_aatk_multithreading 
    PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
               SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
install(TARGETS _aatk_multithreading
        DESTINATION python/pyaatk/multithreading
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/pyaatk/multithreading
)
