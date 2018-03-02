#include <iostream>
#include <cmath>
#include <average-atom-tools/thomas-fermi/atom/potential.h>

int main() {
	AATools::TF::Potential phi;
    phi.setV(1.0);
    phi.setT(1.0);
    phi.setZ(13.0);
    std::cout << phi(1.0) << std::endl;
}