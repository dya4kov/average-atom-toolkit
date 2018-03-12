#pragma once
#include <vector>
#include <cstddef>
#include <average-atom-toolkit/thomas-fermi/quantum-exchange/table.h>

namespace aatk {
namespace TF {
namespace QE {

class Potential {
public:
    Potential();
    Potential(const Potential& potential);
    Potential& operator=(const Potential& potential);
    
    double operator()(const double& x);
    double dx(const double& x);
    
    std::vector<double>& operator()(const std::vector<double>& _xarray);
    std::vector<double>& dx(const std::vector<double>& _xarray);

    double* operator()(const double* x, const std::size_t& n);
    double* dx(const double* x, const std::size_t& n);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setVTZ(
        const double& V, 
        const double& T, 
        const double& Z
    );

    void setTolerance(const double& eps);

private:

    void update();

    double V, T, Z, mu;
    double psi1;
    double tolerance;
    bool   needUpdate;

    const double bestTolerance;
    Table table;
};

}
}
}