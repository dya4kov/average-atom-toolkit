#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

namespace aatk {
namespace TF {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::Half;
using ::numtk::specfunc::FD::ThreeHalf;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSF1 {
    public:
    static const Dimension dim    = 3;
    static const Dimension result = 2;
    RHSF1() : V(1.0), T(1.0), mu(1.0) { 
        rhs = std::bind(&RHSF1::rhsT, this, _1, _2, _3); eval_a(); eval_phi0();
    }
    void set_V (const double& _V)  { V = _V; eval_a(); eval_phi0(); }
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 ?
              std::bind(&RHSF1::rhsT,  this, _1, _2, _3) :
              std::bind(&RHSF1::rhsT0, this, _1, _2, _3) ;
    }
    void set_mu(const double& _mu) { mu = _mu; }

    double param() { return -3.0*std::sqrt(2.0)*V/M_PI/M_PI; }

    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }

    private:

    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = (y[0] + x*x*mu);
        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf(phi/(T*x*x));
            dydx[2] = x*x*x*T*std::sqrt(T)*(phi - phi0)*FDhalf(phi/(T*x*x));
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[2] = 2.0/3.0*(phi - phi0)*std::sqrt(phi)*phi;
        }
    }
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu;
        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        dydx[2] = 2.0/3.0*(phi - phi0)*std::sqrt(phi)*phi;
    }
    void eval_a() { 
        a = std::pow(2.0, 7.0/6.0)
          * std::pow(3.0, 2.0/3.0)
          * std::pow(M_PI, -5.0/3.0)
          * std::pow(V, 2.0/3.0);
    }
    void eval_phi0() {
        phi0 = 1.0/std::pow(3.0*V/(4.0*M_PI), 1.0/3.0);
    }
    double a;
    double V;
    double T;
    double mu;
    double phi0;
    FermiDirac<Half> FDhalf;
};

struct RHSF2 {
    public:
    static const Dimension dim = 3;
    static const Dimension result = 2;
    RHSF2() : V(1.0), T(1.0), mu(1.0) { 
        rhs = std::bind(&RHSF2::rhsT, this, _1, _2, _3); eval_a();
    }
    void set_V (const double& _V)  { V = _V; eval_a(); }
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 ?
              std::bind(&RHSF2::rhsT,  this, _1, _2, _3) :
              std::bind(&RHSF2::rhsT0, this, _1, _2, _3) ;
    }
    void set_mu(const double& _mu) { mu = _mu; }

    double param() { return 4.0*std::sqrt(2.0)*V/M_PI/M_PI; }

    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }
        
    private:

    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = (y[0] + x*x*mu);
        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf(phi/(T*x*x));
            dydx[2] = x*x*x*x*x*T*T*std::sqrt(T)*FD3half(phi/(T*x*x));
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[2] = 2.0/5.0*phi*phi*std::sqrt(phi);
        }
    }
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi = y[0] + x*x*mu;
        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        dydx[2] = 2.0/5.0*phi*phi*std::sqrt(phi);
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
    FermiDirac<Half> FDhalf;
    FermiDirac<ThreeHalf> FD3half;
};

}
}
}