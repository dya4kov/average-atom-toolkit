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

struct RHSPotential {
	static const Dimension dim  = 4;

	RHSPotential() : V(1.0), T(1.0), Z(1.0), mu(4.100577730112) {
        rhs = std::bind(&RHSPotential::rhsT, this, _1, _2, _3); eval_a(); eval_r0();
    }

    ~RHSPotential() { delete drho_sh; }

    void set_V   (const double& _V)   { V = _V; eval_a(); eval_r0(); }
    void set_T   (const double& _T)   { T = _T; T12   = std::sqrt(T); T32   = T*T12;
        rhs = T > 1e-10 
                ? std::bind(&RHSPotential::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSPotential::rhsT0, this, _1, _2, _3);
    }
    void set_Z    (const double& _Z)  { Z = _Z; Z13 = std::pow(Z, 1.0/3.0); Z23 = Z13*Z13; Z43 = Z23*Z23; }
    void set_mu   (const double& _mu) { mu = _mu; }
    void set_eb   (const double& _eb) { eb = _eb; }

    void set_drho(const std::vector<double>& x, const std::vector<double>& y) {
    	drho_sh = new Spline(x, y);
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
        double x2    = x*x;
        double x3    = x*x2;
        double x5    = x2*x3;
        
        double phiMU = y[0] + x2*mu; 
        double phiBE = y[0] + x2*eb; 

        dydx.fill(0.0);

        dydx[0] = 2.0*x*y[1];
        dydx[2] = 2.0*x*y[3];
        if (x > 0) {
            dydx[1] = 2.0*a*x3*T32*FDhalf(phiMU/(T*x2));
        	double drho_tf = a*x*T12*FDmhalf(phiMU/(T*x2));
        	dydx[3] = drho_tf*y[2] + 2.0/x*drho_sh->operator()(x);
            // std::cout << "x = " << x << ", drhoTF = " << drho_tf << ", drhoSH = " << 2.0/x*drho_sh->operator()(x) << std::endl;
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phiMU)*phiMU;
            double drho_tf = 2.0*a*std::sqrt(phiMU);
            dydx[3] = drho_tf*y[2];
            if (phiBE > 0.0) {
                double phi = std::min(phiBE, phiMU);
                dydx[3] -= 4.0/3.0*a*std::sqrt(phi)*phi*std::pow(Z, 4.0/3.0);
            }
        }
        // std::cout << "x = " << x << ", y0 = " << y[0] << ", y1 = " << y[1] << ", y[2] = " << y[2] << ", y3 = " << y[3] << std::endl;
    }
    // rhs at T = 0
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2    = x*x;
        double x3    = x*x2;
        double x5    = x2*x3;
        
        double phiMU = y[0] + x2*mu; 
        double phiBE = y[0] + x2*eb; 

        dydx.fill(0.0);

        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phiMU)*phiMU;

        double drho_tf = 2.0*a*std::sqrt(phiMU);

        dydx[2] = 2.0*x*y[3];
        if (x > 0) {
        	dydx[3] = drho_tf*y[2] + 2.0/x*drho_sh->operator()(x);
        }
        else {
        	dydx[3] = drho_tf*y[2];
            if (phiBE > 0.0) {
                double phi = std::min(phiBE, phiMU);
                dydx[3] -= 4.0/3.0*a*std::sqrt(phi)*phi*std::pow(Z, 4.0/3.0);
            }
        }

    }

    double a;
    double mu;
    double eb;
    double V;
    double T, T12, T32;
    double Z, Z13, Z23, Z43;
    double r0;

    FermiDirac<Half>  FDhalf;
    FermiDirac<MHalf> FDmhalf;
    Spline* drho_sh;

};

}
}
}
}