#pragma once
#include <cmath>
#include <functional>
#include <numeric-tools/ODE/types.h>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

namespace AATools {
namespace TF {
namespace ODE {

using ::numtools::ODE::Array;
using ::numtools::ODE::Dimension;

using ::numtools::specfunc::FermiDirac;
using ::numtools::specfunc::FD::Half;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSAction {
	static const Dimension dim    = 3;
    static const Dimension result = 2;

	RHSAction() : V(1.0), T(1.0), mu(1.0), energyArg(1.0), lambdaArg(1.0) {
        rhs = std::bind(&RHSAction::rhsT, this, _1, _2, _3); eval_a();
    }

    void set_V (const double& _V)  { V = _V; eval_a(); }
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 
                ? std::bind(&RHSAction::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSAction::rhsT0, this, _1, _2, _3);
    }
    void set_mu(const double& _mu) { mu = _mu; }
    void set_l (const double& _l)  { lambdaArg = _l; }
    void set_e (const double& _e)  { energyArg = _e; }
		
	void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }
private:
    void eval_a() {
        a = std::pow(2.0, 7.0/6.0)
          * std::pow(3.0, 2.0/3.0)
          * std::pow(M_PI, -5.0/3.0)
          * std::pow(V, 2.0/3.0);
    }
    // actual rhs 
    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
    // rhs at T > 0
    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu; 
        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf(phi/(T*x*x));
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        }
        double p = energyArg*x*x + y[0] - lambdaArg / (x*x);
        dydx[2] = p > 0.0 ? -std::sqrt(p) : 0.0;
    }
    // rhs at T = 0
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu;

        dydx[0]  = 2.0*x*y[1];
        dydx[1]  = 4.0/3.0*a*std::sqrt(phi)*phi;

        double p = energyArg*x*x + y[0] - lambdaArg / (x*x);

        dydx[2] = p > 0.0 ? -std::sqrt(p) : 0.0;
    }

	double a;
	double energyArg;
	double lambdaArg;
    double mu;
    double T;
    double V;
    FermiDirac<Half> FDhalf;
};

}
}
}