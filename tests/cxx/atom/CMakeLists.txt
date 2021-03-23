add_executable(self-consistent-atom self-consistent-atom.cxx)
target_link_libraries(self-consistent-atom aatk-atom)

install(
	TARGETS self-consistent-atom
	DESTINATION tests/cxx/atom
)