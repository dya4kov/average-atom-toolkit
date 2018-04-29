#pragma once
#include <vector>
#include <cstddef>
#include <functional>

namespace aatk {
namespace TF {
namespace qx {

class ChemicalPotential {
public:
    ChemicalPotential();
    ChemicalPotential(const ChemicalPotential& mu);
    ChemicalPotential& operator=(const ChemicalPotential& mu);

    double operator()(const double& VZ, const double& TZ);
    std::vector<double> operator() (const std::vector<double>& VZ, const std::vector<double>& TZ);
    double* operator()(const double* VZ, const double* TZ, const std::size_t& vsize, const std::size_t& tsize);

    void setZ(const double& Z);
    void setThreadsLimit(const std::size_t& Nthreads);
    void setTolerance(const double& t);

private:

    double tolerance;
    double Z;
    std::size_t threadsLimit;

    void M(const double& V, const double& T, double& result, bool& finished);

    double mu1(const double& V1, const double& T1, const double& M1, const double& tol);

    double mu1_approx(const double& lgV1  /*    T = 0     */);
    double mu1_approx(const double& lgV1, const double& lgT1);

    void updateThreads(
        std::size_t& threads, 
        std::size_t& current, 
        std::size_t& last, 
        bool* finished
    );

    double* evaluate(
        ::std::function<void(const double&, const double&, double&, bool&)> func,
          const double* V, 
          const double* T,
          const std::size_t& vsize, 
          const std::size_t& tsize
    );

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