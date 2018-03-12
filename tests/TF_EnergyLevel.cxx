#include <iostream>
#include <chrono>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

int main() {
    auto start = std::chrono::system_clock::now();

    aatk::TF::EnergyLevel e;
    e.setZ(1.0);
    e.setThreadsLimit(8);
    e.prepareLevelsBelow(15);
        
    for (int n = 1; n < 15; ++n) {
        for (int l = 0; l < n; ++l) {
            std::cout << "e[" << n << "][" << l << "] = " << e[n][l] << std::endl;
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time: " << elapsed.count() << std::endl;
}