#pragma once
#include <vector>
#include <cstddef>
#include <average-atom-tools/thomas-fermi/thermodynamics/chemical-potential.h>

namespace AATools {
namespace TF {

class ElectronDensity {
public:
    ElectronDensity();
    ElectronDensity(const ElectronDensity& eDens);
    ElectronDensity& operator=(const ElectronDensity& eDens);
    
    double operator()(const double& x);
    std::vector<double>& operator()(const std::vector<double>& x);
    double* operator()(const double* x, const std::size_t& n);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setTolerance(const double& eps);

private:

    double V, T, Z;
    double tolerance;

    ChemicalPotential mu;
};

}
}