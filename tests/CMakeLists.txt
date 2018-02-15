add_executable(FDtest FDtest.cxx)
target_link_libraries(FDtest ${NUM_TOOLS_LIBS})
install(TARGETS FDtest DESTINATION tests)

add_executable(Yfunction Yfunction.cxx)
target_link_libraries(Yfunction ${NUM_TOOLS_LIBS})
install(TARGETS Yfunction DESTINATION tests)

add_executable(ODE_PD5test ODE_PD5test.cxx)
install(TARGETS ODE_PD5test DESTINATION tests)

add_executable(ODE_PD853test ODE_PD853test.cxx)
install(TARGETS ODE_PD853test DESTINATION tests)

add_executable(TF_Atom TF_Atom.cxx)
target_link_libraries(TF_Atom ${AVERAGE_ATOM_LIBS})
install(TARGETS TF_Atom DESTINATION tests)

add_executable(TF_Potential TF_Potential.cxx)
target_link_libraries(TF_Potential ${AVERAGE_ATOM_LIBS})
install(TARGETS TF_Potential DESTINATION tests)

add_executable(TF_PotQE TF_PotQE.cxx)
target_link_libraries(TF_PotQE ${AVERAGE_ATOM_LIBS})
install(TARGETS TF_PotQE DESTINATION tests)

add_executable(TF_EnergyLevel TF_EnergyLevel.cxx)
target_link_libraries(TF_EnergyLevel ${AVERAGE_ATOM_LIBS})
install(TARGETS TF_EnergyLevel DESTINATION tests)

add_executable(TF_eStates TF_eStates.cxx)
target_link_libraries(TF_eStates ${AVERAGE_ATOM_LIBS})
install(TARGETS TF_eStates DESTINATION tests)

add_executable(TF_FreeEnergy TF_FreeEnergy.cxx)
target_link_libraries(TF_FreeEnergy ${AVERAGE_ATOM_LIBS})
install(TARGETS TF_FreeEnergy DESTINATION tests)

add_executable(TFQE_FreeEnergy TFQE_FreeEnergy.cxx)
target_link_libraries(TFQE_FreeEnergy ${AVERAGE_ATOM_LIBS})
install(TARGETS TFQE_FreeEnergy DESTINATION tests)