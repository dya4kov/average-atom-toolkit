set(MODULE_NAME "_${PY_PROJECT_NAME}_${MODULE_DIR_NAME}_atom")
# add python module
pybind11_add_module(${MODULE_NAME} ${MODULE_NAME}.cxx)
target_link_libraries(${MODULE_NAME} 
    PRIVATE pybind11::module
            aatk-thomas-fermi
)
install(TARGETS ${MODULE_NAME}
        DESTINATION python/${PY_PROJECT_NAME}/${MODULE_DIR_NAME}/atom
)

add_subdirectory(qx)
add_subdirectory(shell)