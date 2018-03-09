#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

namespace aatk {
namespace TF {
namespace shell {
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

struct RHSdM {

    static const Dimension dim = 3;
    RHSdM() : V(1.0), T(1.0), mu(1.0) {
        rhs = std::bind(&RHSdM::rhsT, this, _1, _2, _3); eval_a();
    }
    void set_V (const double& _V)  { V = _V; eval_a(); } 
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 
                ? std::bind(&RHSdM::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSdM::rhsT0, this, _1, _2, _3);
    }
    void set_mu(const double& _mu) { mu = _mu; }
    
    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }

private:

	std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;

    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu;
        double x2 = x*x;
        double x3 = x*x2;
        double x5 = x2*x3;
        double sqrtT = std::sqrt(T);
        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x3*sqrtT*T*FDhalf(phi/(T*x2));
            dydx[2] = -x5*sqrtT*FDmhalf(phi/(T*x2));
        }
        else {
            dydx[1] = 4.0/3.0*a*sqrt(phi)*phi;
            dydx[2] = 0.0;
        }
    }

    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu;
        double sqrtphi = std::sqrt(phi);
        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*sqrtphi*phi;
        dydx[2] = -x*x*x*x*sqrtphi;
    }

	void eval_a() {
        a = std::pow(2.0, 7.0/6.0)
          * std::pow(3.0, 2.0/3.0)
          * std::pow(M_PI, -5.0/3.0)
          * std::pow(V, 2.0/3.0);
    }
    double a;
    double T;
    double V;
    double mu;
    FermiDirac<Half>  FDhalf;
    FermiDirac<MHalf> FDmhalf;
};

}
}
}
}