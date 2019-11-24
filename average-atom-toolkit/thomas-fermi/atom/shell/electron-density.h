#pragma once
#include <vector>
#include <cstddef>
#include <functional>
#include <average-atom-toolkit/configuration.h>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

namespace aatk {
namespace TF {
namespace shell {

class ElectronDensity {
public:
    ElectronDensity();
    ElectronDensity(const ElectronDensity& eDens);
    ElectronDensity& operator=(const ElectronDensity& eDens);
    
    double operator()(const double x);
    std::vector<double> operator()(const std::vector<double>& x);
    double* operator()(const double* x, const std::size_t n);
    void operator()(const double* x, double* result, const std::size_t n);

    void setV(const double V);
    void setT(const double T);
    void setZ(const double Z);

    void setVTZ(
        const double V, 
        const double T, 
        const double Z
    );

    void setTolerance(const double eps);
    void setBoundary(const double eb);
    void setEnergyLevels(const ::aatk::TF::EnergyLevel& e);
    void setNmax(const std::size_t nmax);

#ifdef ENABLE_MULTITHREADING
    void setThreadsLimit(const std::size_t Nthreads);
#endif    

private:

    double      V, T, Z, mu;
    double      tolerance;
    double      eBoundary;
    std::size_t nmax;
    std::size_t threadsLimit;
    bool        ready;

    ::aatk::TF::EnergyLevel e;

    void accumulate(const double* x, double* rho, const std::size_t size, const std::size_t ithread);
    std::function<void(const double* x, double* rho, const std::size_t size, const std::size_t ithread)> p_accumulate;

};

}
}
}