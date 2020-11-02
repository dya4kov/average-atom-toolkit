add_subdirectory(python)

# add_executable(units units.cxx)
# add_executable(quantity quantity.cxx)
# add_executable(quantity-array quantity-array.cxx)

add_executable(TF_Atom TF_Atom.cxx)
target_link_libraries(TF_Atom aatk-thomas-fermi)

add_executable(TF_chemPot TF_chemPot.cxx)
target_link_libraries(TF_chemPot aatk-thomas-fermi)

add_executable(TF_Potential TF_Potential.cxx)
target_link_libraries(TF_Potential aatk-thomas-fermi)

add_executable(tfqx_potential tfqx_potential.cxx)
target_link_libraries(tfqx_potential aatk-thomas-fermi)

add_executable(TF_EnergyLevel TF_EnergyLevel.cxx)
target_link_libraries(TF_EnergyLevel aatk-thomas-fermi)

add_executable(TF_eStates TF_eStates.cxx)
target_link_libraries(TF_eStates aatk-thomas-fermi)

add_executable(TF_FreeEnergy TF_FreeEnergy.cxx)
target_link_libraries(TF_FreeEnergy aatk-thomas-fermi)

add_executable(TF_FreeEnergy_2nd_derivatives TF_FreeEnergy_2nd_derivatives.cxx)
target_link_libraries(TF_FreeEnergy_2nd_derivatives aatk-thomas-fermi)

add_executable(tfqx_chemPot tfqx_chemPot.cxx)
target_link_libraries(tfqx_chemPot aatk-thomas-fermi)

add_executable(tfqx_FreeEnergy tfqx_FreeEnergy.cxx)
target_link_libraries(tfqx_FreeEnergy aatk-thomas-fermi)

add_executable(TFSH_chemPot TFSH_chemPot.cxx)
target_link_libraries(TFSH_chemPot aatk-thomas-fermi)

add_executable(TFSH_FreeEnergy TFSH_FreeEnergy.cxx)
target_link_libraries(TFSH_FreeEnergy aatk-thomas-fermi)

add_executable(tfsh_wave_function tfsh_wave_function.cxx)
target_link_libraries(tfsh_wave_function aatk-thomas-fermi)

add_executable(tfsh_density tfsh_density.cxx)
target_link_libraries(tfsh_density aatk-thomas-fermi)

install(TARGETS #units
	            #quantity
	            #quantity-array
                TF_Atom
			    TF_chemPot
			    TF_Potential
			    TF_EnergyLevel
			    TF_eStates
			    TF_FreeEnergy
			    TF_FreeEnergy_2nd_derivatives
			    tfqx_potential
			    tfqx_chemPot
			    tfqx_FreeEnergy
			    TFSH_chemPot
			    TFSH_FreeEnergy
			    tfsh_wave_function
			    tfsh_density
	    DESTINATION tests
)