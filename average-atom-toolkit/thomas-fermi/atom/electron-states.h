#pragma once

#include <vector>
#include <cstddef>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

namespace aatk {
namespace TF {

class ElectronStates {
public:

    ElectronStates();

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);
    void setVTZ(
        const double& V, 
        const double& T, 
        const double& Z
    );
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
    std::vector<double> continuous(const std::vector<double>& e);
    std::vector<double> discrete  (const std::vector<double>& e);
    // number of states below each argument energy e (C-style)
    double* continuous(const double* e, const std::size_t& n);
    double* discrete  (const double* e, const std::size_t& n);

    void setNmax(const int& n);
    void setMuShift(const double& mu);
    void setTolerance(const double& eps);
    
#ifdef ENABLE_MULTITHREADING
    void setThreadsLimit(const std::size_t& N);
#endif

    double pseudoDS(const double& e);
    double pseudoCS(const double& e);

    EnergyLevel eLevel();

    double DV();
    double DT();
    double DM();

private:

    std::vector<double> BEroots(const double& eLeft, const double& eRight);
    
    EnergyLevel e;
    double      V, T, Z, mu;
    double      tolerance;
    double      muShift;
    int         nMax;

};

}
}