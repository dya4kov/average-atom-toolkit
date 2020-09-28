add_library(_aatk_semiclassic_atom MODULE _atom.cxx)
target_link_libraries(_aatk_semiclassic_atom 
	PRIVATE pybind11::module 
	        aatk-semiclassic-atom
	        aatk-multithreading
)
set_target_properties(_aatk_semiclassic_atom 
	PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
	           SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
install(TARGETS _aatk_semiclassic_atom
        DESTINATION python/pyaatk/semiclassic
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/pyaatk/semiclassic
)
