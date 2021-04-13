add_executable(atom atom.cxx)
target_link_libraries(atom aatk-thomas-fermi)

add_executable(chemPot chemPot.cxx)
target_link_libraries(chemPot aatk-thomas-fermi)

add_executable(chemPot-qx chemPot-qx.cxx)
target_link_libraries(chemPot-qx aatk-thomas-fermi)

add_executable(chemPot-shell chemPot-shell.cxx)
target_link_libraries(chemPot-shell aatk-thomas-fermi)

add_executable(density-shell density-shell.cxx)
target_link_libraries(density-shell aatk-thomas-fermi)

add_executable(energyLevel energyLevel.cxx)
target_link_libraries(energyLevel aatk-thomas-fermi)

add_executable(eStates eStates.cxx)
target_link_libraries(eStates aatk-thomas-fermi)

add_executable(freeEnergy freeEnergy.cxx)
target_link_libraries(freeEnergy aatk-thomas-fermi)

add_executable(freeEnergy-derivatives freeEnergy-derivatives.cxx)
target_link_libraries(freeEnergy-derivatives aatk-thomas-fermi)

add_executable(freeEnergy-qx freeEnergy-qx.cxx)
target_link_libraries(freeEnergy-qx aatk-thomas-fermi)

add_executable(freeEnergy-shell freeEnergy-shell.cxx)
target_link_libraries(freeEnergy-shell aatk-thomas-fermi)

add_executable(potential potential.cxx)
target_link_libraries(potential aatk-thomas-fermi)

add_executable(potential-qx potential-qx.cxx )
target_link_libraries(potential-qx aatk-thomas-fermi)

add_executable(wave-function wave-function.cxx)
target_link_libraries(wave-function aatk-thomas-fermi)

install(
	TARGETS atom 
	        chemPot 
	        chemPot-qx 
	        chemPot-shell 
	        density-shell 
	        energyLevel 
	        eStates 
	        freeEnergy 
	        freeEnergy-derivatives
	        freeEnergy-qx
	        freeEnergy-shell
	        potential 
	        potential-qx
	        wave-function
	DESTINATION tests/cxx/thomas-fermi
)