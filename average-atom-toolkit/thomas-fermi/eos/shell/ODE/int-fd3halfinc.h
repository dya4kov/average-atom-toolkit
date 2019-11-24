#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>

namespace aatk  {
namespace TF    {
namespace shell {
namespace ODE   {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::Half;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::ThreeHalf;

using ::numtk::specfunc::FermiDiracInc;
using ::numtk::specfunc::FDI::ThreeHalfInc;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHS_int_fd3halfinc {
    static const Dimension dim    = 3;
    static const Dimension result = 2;

    RHS_int_fd3halfinc() : V(1.0), T(1.0), mu(1.0), e(1.0) {
        rhs = std::bind(&RHS_int_fd3halfinc::rhsT, this, _1, _2, _3); eval_a();
    }

    void set_V (const double& _V)  { V = _V; eval_a(); }
    void set_T (const double& _T)  { 
        T   = _T;
        T12 = std::sqrt(T); 
        T32 = T*T12;
        T52 = T*T32; 
        rhs = T > 1e-10 
                ? std::bind(&RHS_int_fd3halfinc::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHS_int_fd3halfinc::rhsT0, this, _1, _2, _3);
    }

    void set_e  (const double& _e  ) { e   = _e;  }
    void set_mu (const double& _mu ) { mu  = _mu; }
    
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
            dydx[1] = 2.0*a*x3*T32*FDhalf(phi/(T*x2));
            if (phiBE <= 0.0) dydx[2] = 0.0;
            else {
                double argX = phi/(T*x2);
                double argY = phiBE/(T*x2);
                if (argY - argX > 25.0)
                    dydx[2] = -x5*T52*FD3half(argX);
                else
                    dydx[2] = -x5*T52*FD3halfI(argX, argY);
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
                dydx[2] = -std::pow(std::min(e, mu)*x2 + y[0], 2.5);
            }
        }
        else {
            dydx[2] = 0.0;
        }
    }

    double a;
    double e;
    double mu;
    double T, T12, T32, T52;
    double V;

    FermiDirac<Half>            FDhalf;
    FermiDirac<ThreeHalf>       FD3half;
    FermiDiracInc<ThreeHalfInc> FD3halfI;
};

}
}
}
}