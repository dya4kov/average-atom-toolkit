#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>

namespace aatk {
namespace TF {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::Half;

using ::numtk::specfunc::FermiDiracInc;
using ::numtk::specfunc::FDI::HalfInc2;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSCS {
    static const Dimension dim    = 3;
    static const Dimension result = 2;

    RHSCS() : V(1.0), T(1.0), mu(1.0), e(1.0), dmu(0.0) {
        rhs = std::bind(&RHSCS::rhsT, this, _1, _2, _3); eval_a();
    }

    void set_V (const double& _V)  { V = _V; eval_a(); }
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 
                ? std::bind(&RHSCS::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSCS::rhsT0, this, _1, _2, _3);
    }

    void set_e  (const double& _e  ) { e   = _e;   }
    void set_mu (const double& _mu ) { mu  = _mu;  }
    void set_dmu(const double& _dmu) { dmu = _dmu; }
    
    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }
private:
    void eval_a() {
        a = std::pow(2.0,  7.0/6.0)
          * std::pow(3.0,  2.0/3.0)
          * std::pow(M_PI,-5.0/3.0)
          * std::pow(V,    2.0/3.0);
    }
    // actual rhs 
    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
    // rhs at T > 0
    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2    = x*x;
        double x3    = x2*x;
        double x5    = x2*x3;
        double phi   = y[0] + x2*mu;
        double phiBE = y[0] + x2*e;
        dydx[0] = 2.0*x*y[1];
        if (x > 0.0) {
            dydx[1] = 2.0*a*x3*std::sqrt(T)*T*FDhalf(phi/(T*x2));
            if (phiBE <= 0.0) dydx[2] = 0.0;
            else {
                double argX = phi/(T*x2) + dmu/T;
                double argY = phiBE/(T*x2);
                if (argY - argX > 25.0)
                    dydx[2] = -x5*FDhalf(argX);
                else
                    dydx[2] = -x5*FDhalfI(argX, argY);
            }
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[2] = 0.0;
        }
    }
    // rhs at T = 0
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2    = x*x;
        double phi   = y[0] + x2*mu;
        double phiBE = y[0] + x2*e;
        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;

        if (x > 0.0) {
            if (phiBE <= 0) dydx[2] = 0.0;
            else {
                dydx[2] = -x2*std::pow((std::min(e, mu) + dmu)*x2 + y[0], 1.5);
            }
        }
        else {
            dydx[2] = 0.0;
        }
    }

    double a;
    double e;
    double mu;
    double dmu;
    double T;
    double V;
    FermiDirac<Half> FDhalf;
    FermiDiracInc<HalfInc2> FDhalfI;
};

struct RHSCSF {
    static const Dimension dim    = 3;
    static const Dimension result = 2;

    RHSCSF() : V(1.0), T(1.0), mu(1.0), dmu(0.0) {
        rhs = std::bind(&RHSCSF::rhsT, this, _1, _2, _3); eval_a();
    }

    void set_V (const double& _V)  { V = _V; eval_a(); }
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 
                ? std::bind(&RHSCSF::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSCSF::rhsT0, this, _1, _2, _3);
    }
    void set_mu(const double& _mu)   {  mu = _mu;  }
    void set_dmu(const double& _dmu) { dmu = _dmu; }
    
    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }
private:
    void eval_a() {
        a = std::pow(2.0,  7.0/6.0)
          * std::pow(3.0,  2.0/3.0)
          * std::pow(M_PI,-5.0/3.0)
          * std::pow(V,    2.0/3.0);
    }
    // actual rhs 
    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
    // rhs at T > 0
    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2  = x*x;
        double x3  = x2*x;
        double x5  = x2*x3;
        double phi = y[0] + x2*mu;

        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x3*std::sqrt(T)*T*FDhalf(phi/(T*x2));
            dydx[2] = -x5*FDhalf(phi/(T*x2) + dmu/T);
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[2] = 0;
        }
    }
    // rhs at T = 0
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2  = x*x;
        double phi = y[0] + x2*mu;
        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        dydx[2] = -x2*std::pow(phi + x2*dmu, 1.5);
    }

    double a;
    double mu;
    double dmu;
    double T;
    double V;
    FermiDirac<Half> FDhalf;
};

}
}
}