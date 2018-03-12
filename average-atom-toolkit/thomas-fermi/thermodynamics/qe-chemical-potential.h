#pragma once
#include <vector>
#include <cstddef>
#include <functional>

namespace aatk {
namespace TF {
namespace QE {

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

};

}
}
}