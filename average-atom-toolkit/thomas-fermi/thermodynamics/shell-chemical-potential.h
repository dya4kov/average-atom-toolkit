#pragma once
#include <vector>
#include <cstddef>

namespace aatk {
namespace TF {
namespace shell {

class ChemicalPotential {
public:
    ChemicalPotential();
    ChemicalPotential(const ChemicalPotential& rps);
    ChemicalPotential& operator=(const ChemicalPotential& rps);

    double operator()(const double& V, const double& T);
    std::vector<double> operator() (const std::vector<double>& VZ, const std::vector<double>& TZ);
    double* operator()(const double* VZ, const double* TZ, const std::size_t& vsize, const std::size_t& tsize);

    void setZ(const double& Z);
    void setThreadsLimit(const std::size_t& N);
    void setTolerance(const double& t);

private:

    double tolerance;
    double Z;
    std::size_t threadsLimit;

    double M(const double& V, const double& T);

};

}
}
}