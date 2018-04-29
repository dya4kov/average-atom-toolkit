#pragma once
#include <cmath>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

namespace aatk {
namespace TF {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::MHalf;
using ::numtk::specfunc::FD::Half;
using ::numtk::specfunc::FD::ThreeHalf;

struct RHSFD2T {
    public:
    static const Dimension dim = 5;
    static const Dimension result = 4;
    RHSFD2T() : V(1.0), T(1.0), mu(1.0) { eval_a(); }
    void set_V (const double& _V)  { V = _V; eval_a(); } 
    void set_T (const double& _T)  { T = _T; }
    void set_mu(const double& _mu) { mu = _mu; }
    double param() { return 3.0*std::sqrt(2.0/T)*V/(T*M_PI*M_PI); }
    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) {
        dydx[0] = 2.0*x*y[1];
        dydx[2] = 2.0*x*y[3];
        double x2  = x*x;
        double x3  = x2*x;
        double x5  = x2*x3;
        double phi = y[0] + x2*mu;
        double sqrtT = std::sqrt(T);
        if (x > 0) {
            dydx[1] = 2.0*a*x3*sqrtT*T*FDhalf(phi/(T*x2));
            dydx[3] = a*x*sqrtT*(
                               3.0*x2*FDhalf(phi/(T*x2))
                    + (y[2] - phi/T)*FDMhalf(phi/(T*x2))
                    );
            dydx[4] = 5.0*x5*T*T*FD3half(phi/(T*x2)) 
                    + 3.0*x3*T*(y[2]*T - 2.0*phi)*FDhalf(phi/(T*x2))
                    + x*phi*(phi - y[2]*T)*FDMhalf(phi/(T*x2));
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[3] = 2.0*a*std::sqrt(phi)*y[2];
            dydx[4] = 0.0;
        }
    }
    private:
    void eval_a() {
        a = std::pow(2.0, 7.0/6.0)
          * std::pow(3.0, 2.0/3.0)
          * std::pow(M_PI, -5.0/3.0)
          * std::pow(V, 2.0/3.0);
    }
    double a;
    double V;
    double T;
    double mu;
    FermiDirac<MHalf>     FDMhalf;
    FermiDirac<Half>      FDhalf;
    FermiDirac<ThreeHalf> FD3half;
};

}
}
}