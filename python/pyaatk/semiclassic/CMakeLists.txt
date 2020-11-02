# use folder name as library name
get_filename_component(MODULE_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(MODULE_NAME "_${PY_PROJECT_NAME}_${MODULE_DIR_NAME}")
# add python module
pybind11_add_module(${MODULE_NAME} ${MODULE_NAME}.cxx)
target_link_libraries(${MODULE_NAME} 
    PRIVATE pybind11::module 
            aatk-semiclassic
)
install(TARGETS ${MODULE_NAME}
        DESTINATION python/${PY_PROJECT_NAME}/${MODULE_DIR_NAME}
)
