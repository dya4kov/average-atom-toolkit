#pragma once

#include <vector>
#include <cstddef>

#include <average-atom-tools/thomas-fermi/atom/potential.h>
#include <average-atom-tools/thomas-fermi/atom/energy-level.h>

namespace AATools {
namespace TF {

class ElectronStates {
public:

    ElectronStates();

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);
    // fill numbers for discrete levels
    double operator()(const int& n);
    double operator()(const int& n, const int& l);
    // boundary energy between discrete and continuous states
    double eBoundary();
    // total number of states
    double continuous();
    double discrete();
    // number of states below argument energy e
    double continuous(const double& e);
    double discrete  (const double& e);
    // number of states below each argument energy e
    std::vector<double>& continuous(const std::vector<double>& e);
    std::vector<double>& discrete  (const std::vector<double>& e);
    // number of states below each argument energy e (C-style)
    double* continuous(const double* e, const std::size_t& n);
    double* discrete  (const double* e, const std::size_t& n);

    void setNmax(const int& n);
    void setMuShift(const double& mu);
    void setTolerance(const double& eps);

    double pseudoDS(const double& e);
    double pseudoCS(const double& e);

private:

    std::vector<double> BEroots(const double& eLeft, const double& eRight);
    
    EnergyLevel e;
    Potential   phi;

    double      V1, T1, mu1;
    double      VZ, TZ, muZ;
    double      tolerance;
    double      muShift;

    int         nMax;

};

}
}