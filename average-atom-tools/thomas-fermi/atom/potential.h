#pragma once
#include <vector>
#include <cstddef>
#include <average-atom-tools/thomas-fermi/thermodynamics/chemical-potential.h>

namespace AATools {
namespace TF {

class Potential {
public:
    Potential();
    Potential(const Potential& potential);
    Potential& operator=(const Potential& potential);
    
    // double mu();

    double operator()(const double& x);
    double dx(const double& x);
    
    std::vector<double>& operator()(const std::vector<double>& x);
    std::vector<double>& dx(const std::vector<double>& x);

    double* operator()(const double* x, const std::size_t& n);
    double* dx(const double* x, const std::size_t& n);

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