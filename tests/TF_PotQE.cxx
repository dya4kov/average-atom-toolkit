#include <iostream>
#include <cmath>
#include <average-atom-tools/thomas-fermi/quantum-exchange/potential.h>

int main() {
	AATools::TF::QE::Potential psi;
    psi.setV(1.32455);
    psi.setT(1.32455);
    std::cout << psi(1.0) << std::endl;
}