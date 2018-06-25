#include <iostream>
#include <chrono>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

int main() {
    auto start = std::chrono::system_clock::now();

    aatk::TF::EnergyLevel e;
    e.setVTZ(100.0, 1.0, 47.0);
    // e.prepareLevelsBelow(15);

    e.setTolerance(1e-7);

    for (int n = 1; n < 15; ++n) {
        for (int l = 0; l < n; ++l) {
            std::cout << "e[" << n << "][" << l << "] = " << e(n,l) << std::endl;
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time: " << elapsed.count() << std::endl;
}