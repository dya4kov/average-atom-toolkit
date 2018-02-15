#include <iostream>
#include <cmath>
#include <chrono>
#include <average-atom-tools/thomas-fermi/thermodynamics/qe-free-energy.h>

int main() {
    AATools::TF::QE::FreeEnergy F;
    F.setThreadsLimit(6);
    int nV = 5;
    int nT = 5;
    // double* V = new double[nV];
    // double* T = new double[nT];
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

    // auto result = F.D2T(V, T, nV, nT);
    auto result = F.DT(V, T);

    for (int iV = 0; iV < nV; ++iV) {
        for (int iT = 0; iT < nT; ++iT) {
            std::cout << result[iV*nT + iT] << " ";
        }
        std::cout << std::endl;
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time parallel: " << elapsed.count() << std::endl;

    //start = std::chrono::system_clock::now();

    //for (int iV = 0; iV < nV; ++iV) {
    //    for (int iT = 0; iT < nT; ++iT) {
    //        F(V[iV], T[iT]);
    //    }
    //}

    //end = std::chrono::system_clock::now();
    //elapsed = end - start;
    //std::cout << "elapsed time serial: " << elapsed.count() << std::endl;

}