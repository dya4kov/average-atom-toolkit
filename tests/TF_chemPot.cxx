#include <iostream>
#include <cmath>
#include <average-atom-tools/thomas-fermi/thermodynamics/chemical-potential.h>

int main() {
    AATools::TF::ChemicalPotential mu;
    mu.setZ(1.0);
    int nV = 11;
    int nT = 11;
    std::vector<double> V(nV);
    std::vector<double> T(nT);

    double Vmin = 1.2341; double Vmax = 999.1231;
    double Tmin = 0.1; double Tmax = 999.1231;

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
    std::cout << std::endl;
    for (int iV = 0; iV < nV; ++iV) {
        for (int iT = 0; iT < nT; ++iT) {
            std::cout << mu(V[iV], T[iT]) << " ";
        }
        std::cout << std::endl;
    }

}