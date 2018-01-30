add_executable(FDtest FDtest.cxx)
target_link_libraries(FDtest ${NUM_TOOLS_LIBS})
install(TARGETS FDtest DESTINATION tests)