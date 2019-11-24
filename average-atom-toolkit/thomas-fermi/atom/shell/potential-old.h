#pragma once
#include <vector>
#include <cstddef>

#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>

namespace aatk {
namespace TF {
namespace shell {

class Potential {
public:
    Potential();
    Potential(const Potential& potential);
    Potential& operator=(const Potential& potential);
    
    double operator()(const double& x);
    double dx(const double& x);
    
    std::vector<double> operator()(const std::vector<double>& _xarray);
    std::vector<double> dx(const std::vector<double>& _xarray);

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
    void setNmax(const int& nmax);

private:

    void prepare();

    double V, T, Z, mu, dmush;
    double tolerance;
    double eBoundary;
    bool   ready;
    int    nmax;

    ::aatk::TF::ElectronStates N;

    std::vector<double> e;
    std::vector<double> rpi;
    std::vector<double> rpo;
    std::vector<double> ksi0;
    std::vector<double> ksi21;
    std::vector<double> sign;
};

}
}
}