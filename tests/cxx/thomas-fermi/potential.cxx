#include <iostream>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom/potential.h>

int main() {
	aatk::TF::Potential phi;
    phi.setV(1.0);
    phi.setT(1.0);
    phi.setZ(3.0);
    std::cout << phi(0.0) << std::endl;
}