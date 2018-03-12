add_subdirectory(atom)
add_subdirectory(thermodynamics)
install(FILES Atom.py DESTINATION python/pyaatk/TF)
install(FILES __init__.py DESTINATION python/pyaatk/TF)