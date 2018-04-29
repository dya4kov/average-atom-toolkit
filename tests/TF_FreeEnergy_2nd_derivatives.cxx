#include <iostream>
#include <cmath>
#include <chrono>
#include <average-atom-toolkit/thomas-fermi/eos/free-energy.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

int main() {
	aatk::TF::FreeEnergy        F;
	aatk::TF::ChemicalPotential M;

	double T = 1.0;
	double V = 1.0;
	double Z = 1.0;

	double dV = 0.01;
	double dT = 0.01;

	M.setZ(Z);
	F.setZ(Z);
	F.setTolerance(1e-9);

	for (double T = 1.0; T < 10.1; T += 1.0) {
		// std::cout << (F.DT(V, T + dT) - F.DT(V, T - dT))/(2.0*dT) << "   ";
		std::cout << F.D2T(V, T) << "   ";
		std::cout << F.D2T2(V, T) << std::endl;
	}

	return 0;
}