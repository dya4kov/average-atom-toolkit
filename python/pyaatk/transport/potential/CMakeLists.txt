add_library(_aatk_transport_potential MODULE _potential.cxx)
target_link_libraries(_aatk_transport_potential 
	PRIVATE pybind11::module 
	        aatk-transport 
)
set_target_properties(_aatk_transport_potential 
	PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
	           SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
install(TARGETS _aatk_transport_potential
        DESTINATION python/pyaatk/transport/potential
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/pyaatk/transport/potential
)

