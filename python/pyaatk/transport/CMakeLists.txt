add_subdirectory(potential)
add_library(_aatk_transport_cross_section MODULE _cross-section.cxx)
target_link_libraries(_aatk_transport_cross_section 
	PRIVATE pybind11::module 
	        aatk-transport 
)
set_target_properties(_aatk_transport_cross_section 
	PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
	           SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
install(TARGETS _aatk_transport_cross_section
        DESTINATION python/pyaatk/transport
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/pyaatk/transport
)

