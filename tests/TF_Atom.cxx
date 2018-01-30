#include <iostream>
#include <average-atom-tools/thomas-fermi/atom.h>

int main() {
	AATools::TF::Atom a;
    double V = 1.0;
    double T = 1.0;
    double Z = 1.0;
    a.setVTZ(V, T, Z);
    std::cout << a.e(1,0) << std::endl;
    std::cout << a.e(2,0) << std::endl;
    std::cout << a.e(2,1) << std::endl;
    std::cout << a.e(3,0) << std::endl;
    std::cout << a.e(3,1) << std::endl;
    std::cout << a.e(3,2) << std::endl;
    std::cout << a.e(4,0) << std::endl;
    std::cout << a.e(4,1) << std::endl;
    std::cout << a.e(4,2) << std::endl;
    std::cout << a.e(4,3) << std::endl;
}