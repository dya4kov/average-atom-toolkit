#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <average-atom-toolkit/thomas-fermi/eos/shell/chemical-potential.h>

int main() {
    aatk::TF::shell::ChemicalPotential mu;
    mu.setZ(1.0);
    int nV = 6;
    int nT = 6;
    std::vector<double> V(nV);
    std::vector<double> T(nT);

    double Vmin = 1.0; double Vmax = 10.0;
    double Tmin = 1.0; double Tmax = 10.0;

    for (int iV = 0; iV < nV; ++iV)  {
        V[iV] = Vmin + iV*(Vmax - Vmin)/(nV - 1);
    }
    for (int iT = 0; iT < nT; ++iT) {
        T[iT] = Tmin + iT*(Tmax - Tmin)/(nT - 1);
    }

    auto start = std::chrono::system_clock::now();

    auto result = mu(V, T);

    for (int iV = 0; iV < nV; ++iV) {
        for (int iT = 0; iT < nT; ++iT) {
            std::cout << result[iV*nT + iT] << " ";
            // std::cout << mu(V[iV], T[iT])   << " ";
        }
        std::cout << std::endl;
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time parallel: " << elapsed.count() << std::endl;
}