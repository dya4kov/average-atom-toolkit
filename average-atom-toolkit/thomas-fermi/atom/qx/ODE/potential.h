#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/Yfunction.h>

namespace aatk {
namespace TF {
namespace qx {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::Half;
using ::numtk::specfunc::FD::MHalf;
using ::numtk::specfunc::Yfunction;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSPotential {
    public:
    static const Dimension dim = 4;
    RHSPotential() : V(1.0), T(1.0), mu(1.0) { 
        rhs = std::bind(&RHSPotential::rhsT, this, _1, _2, _3); eval_a();
    }
    void set_V (const double _V)  { V = _V; eval_a(); }
    void set_T (const double _T)  { T = _T; 
        rhs = T > 1e-10 
                ? std::bind(&RHSPotential::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSPotential::rhsT0, this, _1, _2, _3);
    }
    void set_mu(const double _mu) { mu = _mu; }
    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }

    private:

    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        dydx[0] = 2.0*x*y[1];
        dydx[2] = 2.0*x*y[3];
        if (x > 0) {
            dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf((y[0] + x*x*mu)/(T*x*x));
            dydx[3] = 2*x*a*std::sqrt(T)*(0.5*FDmhalf((y[0] + x*x*mu)/(T*x*x))*y[2] 
                    +   std::sqrt(T)*x*x*Y.derivative((y[0] + x*x*mu)/(T*x*x)));
        }
        else {
            double phi = y[0] + x*x*mu;
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[3] = 2.0*a*(std::sqrt(phi)*y[2] + 22.0/3.0*x*phi);
        }
    }
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu;
        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        dydx[2] = 2.0*x*y[3];
        dydx[3] = 2.0*a*(std::sqrt(phi)*y[2] + 22.0/3.0*x*phi);
    }

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
    FermiDirac<Half>  FDhalf;
    FermiDirac<MHalf> FDmhalf;
    Yfunction Y;
};
}
}
}
}