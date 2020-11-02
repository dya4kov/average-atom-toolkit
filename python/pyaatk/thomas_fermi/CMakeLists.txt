# use folder name as library name
get_filename_component(MODULE_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_subdirectory(atom)
add_subdirectory(eos)

install(FILES Atom.py DESTINATION python/${PY_PROJECT_NAME}/${MODULE_DIR_NAME})