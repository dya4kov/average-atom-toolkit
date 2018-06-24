#pragma once
#include <vector>
#include <cstddef>

#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>

namespace aatk {
namespace TF {
namespace shell {

class WaveFunction {
public:
    WaveFunction();
    WaveFunction(const WaveFunction& wf);
    WaveFunction& operator=(const WaveFunction& wf);

    double norm(const double& e, const double& lambda);
    
    double operator()(const double& e, const double& lambda, const double& x);
    std::vector<double> operator()(const double& e, const double& lambda, const std::vector<double>& x);
    double* operator()(const double& e, const double& lambda, const double* x, const std::size_t& size);

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

    void setNorm();
    void setParam(const double& e, const double& l);

    double V, T, Z, mu;
    double tolerance;
    double energy, lambda;

    double normValue; 
    bool   ready;

    RotatePoints RP;
};

}
}
}