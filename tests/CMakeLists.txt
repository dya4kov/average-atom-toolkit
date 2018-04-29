add_executable(units units.cxx)
add_executable(quantity quantity.cxx)
add_executable(quantity-array quantity-array.cxx)

add_executable(FDtest FDtest.cxx)
target_link_libraries(FDtest ${NUM_TOOLS_LIBS})

add_executable(bessel bessel_test.cxx)
target_link_libraries(bessel ${NUM_TOOLS_LIBS})

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

add_executable(TF_EnergyLevel TF_EnergyLevel.cxx)
target_link_libraries(TF_EnergyLevel ${AVERAGE_ATOM_LIBS})

add_executable(TF_eStates TF_eStates.cxx)
target_link_libraries(TF_eStates ${AVERAGE_ATOM_LIBS})

add_executable(TF_FreeEnergy TF_FreeEnergy.cxx)
target_link_libraries(TF_FreeEnergy ${AVERAGE_ATOM_LIBS})

add_executable(TF_FreeEnergy_2nd_derivatives TF_FreeEnergy_2nd_derivatives.cxx)
target_link_libraries(TF_FreeEnergy_2nd_derivatives ${AVERAGE_ATOM_LIBS})

add_executable(tfqx_chemPot tfqx_chemPot.cxx)
target_link_libraries(tfqx_chemPot ${AVERAGE_ATOM_LIBS})

add_executable(tfqx_FreeEnergy tfqx_FreeEnergy.cxx)
target_link_libraries(tfqx_FreeEnergy ${AVERAGE_ATOM_LIBS})

add_executable(TFSH_chemPot TFSH_chemPot.cxx)
target_link_libraries(TFSH_chemPot ${AVERAGE_ATOM_LIBS})

add_executable(TFSH_FreeEnergy TFSH_FreeEnergy.cxx)
target_link_libraries(TFSH_FreeEnergy ${AVERAGE_ATOM_LIBS})

add_executable(tfsh_wave_function tfsh_wave_function.cxx)
target_link_libraries(tfsh_wave_function ${AVERAGE_ATOM_LIBS})

add_executable(tfsh_density tfsh_density.cxx)
target_link_libraries(tfsh_density ${AVERAGE_ATOM_LIBS})

install(TARGETS units
	            quantity
	            quantity-array
                FDtest
                bessel
			    Yfunction
			    ODE_PD5test
			    ODE_PD853test
			    TF_Atom
			    TF_chemPot
			    TF_Potential
			    TF_EnergyLevel
			    TF_eStates
			    TF_FreeEnergy
			    TF_FreeEnergy_2nd_derivatives
			    tfqx_chemPot
			    tfqx_FreeEnergy
			    TFSH_chemPot
			    TFSH_FreeEnergy
			    tfsh_wave_function
			    tfsh_density
	    DESTINATION tests
)