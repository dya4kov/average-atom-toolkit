#pragma once
#include <cmath>
#include <functional>
#include <numeric-tools/ODE/types.h>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

namespace AATools {
namespace TF {
namespace shell {
namespace ODE {

using ::numtools::ODE::Array;
using ::numtools::ODE::Dimension;

using ::numtools::specfunc::FermiDirac;
using ::numtools::specfunc::FD::Half;
using ::numtools::specfunc::FD::MHalf;
using ::numtools::specfunc::Yfunction;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSEnergyInt {

    static const Dimension dim = 3;
    RHSEnergyInt() : V(1.0), T(1.0), mu(1.0) {
        rhs = std::bind(&RHSEnergyInt::rhsT, this, _1, _2, _3); eval_a();
    }
    void set_V (const double& _V)  { V = _V; eval_a(); }
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 
                ? std::bind(&RHSEnergyInt::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSEnergyInt::rhsT0, this, _1, _2, _3);
    }
    void set_mu(const double& _mu) { mu = _mu; }

    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }

private:

    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;

    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2 = x*x;
        double x3 = x2*x;
        double phi = y[0] + x2*mu;
        double sqrtT = std::sqrt(T);
        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x3*sqrtT*T*FDhalf(phi/(T*x2));
            dydx[2] = -x3*sqrtT*phi*FDmhalf(phi/(T*x2));
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[2] = 0.0;
        }
    }

    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu;
        double sqrtphi = std::sqrt(phi);
        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*sqrtphi*phi;
        dydx[2] = -x*x*sqrtphi*phi;
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
};
