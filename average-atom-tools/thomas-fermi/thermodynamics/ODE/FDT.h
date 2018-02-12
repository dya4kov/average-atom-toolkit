#pragma once
#include <cmath>
#include <numeric-tools/ODE/types.h>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

namespace AATools {
namespace TF {
namespace ODE {

using ::numtools::ODE::Array;
using ::numtools::ODE::Dimension;

using ::numtools::specfunc::FermiDirac;
using ::numtools::specfunc::FD::Half;
using ::numtools::specfunc::FD::ThreeHalf;

struct RHSFDT {
    public:
    static const Dimension dim = 3;
    static const Dimension result = 2;
    RHSFDT() : V(1.0), T(1.0), mu(1.0) { eval_a(); }
    void set_V (const double& _V)  { V = _V; eval_a(); } 
    void set_T (const double& _T)  { T = _T; }
    void set_mu(const double& _mu) { mu = _mu; }
    double param() { return -2.0*std::sqrt(2.0*T)*V*T/(M_PI*M_PI); }
    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) {
        double phi = (y[0] + x*x*mu);
        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x*x*x*T*std::sqrt(T)*FDhalf(phi/(T*x*x));
            dydx[2] = x*x*x*x*x*(3.0*phi/(T*x*x)*FDhalf(phi/(T*x*x)) - 5.0*FD3half(phi/(T*x*x)));
        }
        else { 
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[2] = 0.0;
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
    double b;
    double V;
    double T;
    double mu;
    FermiDirac<Half> FDhalf;
    FermiDirac<ThreeHalf> FD3half;
};

}
}
}