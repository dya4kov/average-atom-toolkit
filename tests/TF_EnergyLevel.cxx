#include <iostream>
#include <chrono>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

int main() {
    aatk::TF::EnergyLevel e;

    e.setZ(1.0);
        
    auto start = std::chrono::system_clock::now();
    for (int n = 1; n < 15; ++n) {
        for (int l = 0; l < n; ++l) {
            std::cout << "e[" << n << "][" << l << "] = " << e[n][l] << std::endl;
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time: " << elapsed.count() << std::endl;
}