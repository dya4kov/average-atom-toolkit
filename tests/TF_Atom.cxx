#include <iostream>
#include <chrono>
#include <cmath>

#include <average-atom-toolkit/thomas-fermi/atom.h>

int main() {
	aatk::TF::Atom atom;
    atom.setZ(13.0);
    auto start = std::chrono::system_clock::now();
    for (int n = 1; n < 50; ++n) {
        for (int l = 0; l < n; ++l) {
            std::cout << "e[" << n << "][" << l << "] = " << atom.e[n][l]*std::pow(13.0, -4.0/3.0) << std::endl;
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time: " << elapsed.count() << std::endl;
}