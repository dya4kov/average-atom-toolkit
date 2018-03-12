add_executable(FDtest FDtest.cxx)
target_link_libraries(FDtest ${NUM_TOOLS_LIBS})

add_executable(Yfunction Yfunction.cxx)
target_link_libraries(Yfunction ${NUM_TOOLS_LIBS})

add_executable(ODE_PD5test ODE_PD5test.cxx)
add_executable(ODE_PD853test ODE_PD853test.cxx)

add_executable(TF_Atom TF_Atom.cxx)
target_link_libraries(TF_Atom ${AVERAGE_ATOM_LIBS})

add_executable(TF_chemPot TF_chemPot.cxx)
target_link_libraries(TF_chemPot ${AVERAGE_ATOM_LIBS})

add_executable(TF_Potential TF_Potential.cxx)
target_link_libraries(TF_Potential ${AVERAGE_ATOM_LIBS})

add_executable(TF_PotQE TF_PotQE.cxx)
target_link_libraries(TF_PotQE ${AVERAGE_ATOM_LIBS})

add_executable(TF_EnergyLevel TF_EnergyLevel.cxx)
target_link_libraries(TF_EnergyLevel ${AVERAGE_ATOM_LIBS})

add_executable(TF_eStates TF_eStates.cxx)
target_link_libraries(TF_eStates ${AVERAGE_ATOM_LIBS})

add_executable(TF_FreeEnergy TF_FreeEnergy.cxx)
target_link_libraries(TF_FreeEnergy ${AVERAGE_ATOM_LIBS})

add_executable(TFQE_chemPot TFQE_chemPot.cxx)
target_link_libraries(TFQE_chemPot ${AVERAGE_ATOM_LIBS})

add_executable(TFQE_FreeEnergy TFQE_FreeEnergy.cxx)
target_link_libraries(TFQE_FreeEnergy ${AVERAGE_ATOM_LIBS})

add_executable(TFSH_chemPot TFSH_chemPot.cxx)
target_link_libraries(TFSH_chemPot ${AVERAGE_ATOM_LIBS})

add_executable(TFSH_FreeEnergy TFSH_FreeEnergy.cxx)
target_link_libraries(TFSH_FreeEnergy ${AVERAGE_ATOM_LIBS})

install(TARGETS FDtest
			    Yfunction
			    ODE_PD5test
			    ODE_PD853test
			    TF_Atom
			    TF_chemPot
			    TF_Potential
			    TF_PotQE
			    TF_EnergyLevel
			    TF_eStates
			    TF_FreeEnergy
			    TFQE_chemPot
			    TFQE_FreeEnergy
			    TFSH_chemPot
			    TFSH_FreeEnergy
	    DESTINATION tests
)