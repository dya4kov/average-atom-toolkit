#include <iostream>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/chemical-potential.h>

int main() {
    aatk::TF::ChemicalPotential mu;
    mu.setZ(1.0);
    mu.setThreadsLimit(2);
    int nV = 11;
    int nT = 11;
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

    // std::cout << mu(1.0, 2.0) << std::endl;

    auto result = mu(V, T);

    for (int iV = 0; iV < nV; ++iV) {
        for (int iT = 0; iT < nT; ++iT) {
            std::cout << result[iV*nT + iT] << " ";
        }
        std::cout << std::endl;
    }
}