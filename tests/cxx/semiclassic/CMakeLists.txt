add_executable(self-consistent self-consistent.cxx)
target_link_libraries(self-consistent aatk-semiclassic)
install(
	TARGETS self-consistent
	DESTINATION tests/cxx/semiclassic
)