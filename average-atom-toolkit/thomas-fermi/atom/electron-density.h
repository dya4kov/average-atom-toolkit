#pragma once
#include <vector>
#include <cstddef>

namespace aatk {
namespace TF {

class ElectronDensity {
public:
    ElectronDensity();
    // ElectronDensity(const ElectronDensity& eDens);
    // ElectronDensity& operator=(const ElectronDensity& eDens);
    
    double operator()(const double& x);
    std::vector<double> operator()(const std::vector<double>& x);
    double* operator()(const double* x, const std::size_t& n);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setVTZ(
        const double& V, 
        const double& T, 
        const double& Z
    );

    void setTolerance(const double& eps);
    void setBoundary(const double& eb);

private:

    double V, T, Z, mu;
    double tolerance;
    double eb;
};

}
}