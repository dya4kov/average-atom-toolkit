#include <iostream>
#include <chrono>
#include <cmath>
#include <average-atom-tools/thomas-fermi/atom/electron-states.h>

int main() {
    AATools::TF::ElectronStates N;
    N.setV(200.0);
    N.setT(50.0);
    N.setZ(13.0);
    auto start = std::chrono::system_clock::now();
    std::cout << N.continuous() << std::endl;
    std::cout << N.discrete() << std::endl;
    double eBoundary = N.eBoundary();
    std::cout << eBoundary << std::endl;
    std::cout << N.continuous(eBoundary) << std::endl;
    std::cout << N.discrete(eBoundary) << std::endl;
    for (int n = 1; n < 5; ++n) {
        for (int l = 0; l < n; ++l) {
            std::cout << "N[" << n << "][" << l << "] = " << N(n,l) << std::endl;
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time: " << elapsed.count() << std::endl;
}