#include <iostream>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom/qx/potential.h>

int main() {
	aatk::TF::qx::Potential phi;
    phi.setV(1.0);
    phi.setT(1.0);
    phi.setZ(3.0);
    std::cout << phi(0.0) << std::endl;
}