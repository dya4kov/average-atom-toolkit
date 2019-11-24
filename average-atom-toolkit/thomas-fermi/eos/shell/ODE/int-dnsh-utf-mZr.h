#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/interpolation/spline.h>

namespace aatk {
namespace TF {
namespace shell {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::Half;
using ::numtk::specfunc::FD::MHalf;
using ::numtk::interpolation::Spline;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHS_int_dnsh_utf_mZr {
	static const Dimension dim  = 3;

	RHS_int_dnsh_utf_mZr() : V(1.0), T(1.0), Z(1.0), mu(4.100577730112) {
        rhs = std::bind(&RHS_int_dnsh_utf_mZr::rhsT, this, _1, _2, _3); eval_a(); eval_r0();
    }

    ~RHS_int_dnsh_utf_mZr() { delete dnsh; }

    void set_V   (const double& _V)   { V = _V; eval_a(); eval_r0(); }
    void set_T   (const double& _T)   { 
        T   = _T; 
        T12 = std::sqrt(T); 
        T32 = T*T12;
        T52 = T*T32;
        rhs = T > 1e-10 
                ? std::bind(&RHS_int_dnsh_utf_mZr::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHS_int_dnsh_utf_mZr::rhsT0, this, _1, _2, _3);
    }
    void set_Z   (const double& _Z)  { Z = _Z; Z13 = std::pow(Z, 1.0/3.0); Z23 = Z13*Z13; Z43 = Z23*Z23; }
    void set_mu  (const double& _mu) { mu = _mu; }

    void set_dnsh(const std::vector<double>& x, const std::vector<double>& y) {
    	dnsh = new Spline(x, y);
    }

    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }

private:
    void eval_a() {
        a = std::pow(2.0, 7.0/6.0)
          * std::pow(3.0, 2.0/3.0)
          * std::pow(M_PI, -5.0/3.0)
          * std::pow(V, 2.0/3.0);
    }
    void eval_r0() {
        r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
    }
    // actual rhs 
    std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
    // rhs at T > 0
    void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2  = x*x;
        double x3  = x*x2;
        double x5  = x2*x3;
        
        double phi = y[0] + x2*mu;

        dydx[0] = 2.0*x*y[1];

        if (x > 0) {
            dydx[1] = 2.0*a*x3*T32*FDhalf(phi/(T*x2));
            dydx[2] = -dnsh->operator()(x)*(y[0] - 1.0/r0)/x;
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            dydx[2] = 0.0;
        }

        // std::cout << "x = " << x << ", dnsh = " << dnsh->operator()(x) << std::endl;

    }
    // rhs at T = 0
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2  = x*x;
        double x3  = x*x2;
        double x5  = x2*x3;
        
        double phi = y[0] + x2*mu;

        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        dydx[2] = -dnsh->operator()(x)*(y[0] - 1.0/r0)/x;

    }

    double a;
    double mu;
    double V;
    double T, T12, T32, T52;
    double Z, Z13, Z23, Z43;
    double r0;

    FermiDirac<Half>  FDhalf;
    FermiDirac<MHalf> FDmhalf;
    Spline* dnsh;

};

}
}
}
}