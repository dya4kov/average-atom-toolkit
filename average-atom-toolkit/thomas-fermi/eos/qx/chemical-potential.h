#pragma once
#include <vector>
#include <cstddef>
#include <functional>

#include <average-atom-toolkit/configuration.h>

namespace aatk {
namespace TF {
namespace qx {

class ChemicalPotential {
public:
    ChemicalPotential();
    ChemicalPotential(const ChemicalPotential& mu);
    ChemicalPotential& operator=(const ChemicalPotential& mu);

    double operator()(const double VZ, const double TZ);
    std::vector<double> operator() (const std::vector<double>& VZ, const std::vector<double>& TZ);
    double* operator()(const double* VZ, const double* TZ, const std::size_t vsize, const std::size_t tsize);

    void setZ(const double Z);
    void setTolerance(const double t);

#ifdef ENABLE_MULTITHREADING
    void setThreadsLimit(const std::size_t Nthreads);
#endif

private:

    double tolerance;
    double Z;

    double M(const double V, const double T);

    double mu1(const double V1, const double T1, const double M1, const double tol);

    double mu1_approx(const double lgV1  /*    T = 0    */);
    double mu1_approx(const double lgV1, const double lgT1);

#ifdef ENABLE_MULTITHREADING
    std::size_t threadsLimit;
#endif

    void evaluate(
        ::std::function<double(const double, const double)> func,
          const double* V, 
          const double* T,
                double* result,
          const std::size_t vsize, 
          const std::size_t tsize,
          const std::size_t ithread
    );

    std::function<void(
        ::std::function<double(const double, const double)> func,
          const double* V, 
          const double* T, 
                double* result, 
          const std::size_t vsize, 
          const std::size_t tsize,
          const std::size_t ithread)> 
    p_evaluate;

    // basic buffer
    static std::vector<double>   table;
    static const int    vSize,   tSize;
    static const double lgV0,    lgT0;
    static const double lgVstep, lgTstep;
    static const double bestTolerance;

};

}
}
}