#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/bessel/Jnu.h>
#include <numeric-toolkit/specfunc/bessel/Ynu.h>
#include <numeric-toolkit/specfunc/bessel/Knu.h>

namespace aatk {
namespace TF {
namespace shell {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::Half;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSnorm {
	static const Dimension dim    = 4;
    static const Dimension result = 3;

	RHSnorm() : V(1.0), T(1.0), mu(1.0), energyArg(1.0), lambdaArg(1.0) {
        rhs = std::bind(&RHSnorm::rhsT, this, _1, _2, _3); eval_a(); eval_r0();
    }

    void set_V (const double& _V)  { V = _V; eval_a(); eval_r0(); }
    void set_T (const double& _T)  { T = _T; 
        rhs = T > 1e-10 
                ? std::bind(&RHSnorm::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSnorm::rhsT0, this, _1, _2, _3);
    }
    void set_Z (const double& _Z)  { Z = _Z; Z13 = std::pow(Z, 1.0/3.0); Z23 = Z13*Z13; }
    void set_mu(const double& _mu) { mu = _mu; }
    void set_l (const double& _l)  { lambdaArg = _l; }
    void set_e (const double& _e)  { energyArg = _e; }
    
    void set_ksi0 (const double& _ksi0)  { ksi0  = _ksi0;  }
    void set_ksi21(const double& _ksi21) { ksi21 = _ksi21; }

    void set_RP(const double& _xmin, const double& _xmax) { xmin = _xmin; xmax = _xmax; } 

    void set_sign2(const double& _sign) { sign_2 = _sign; }

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
        double phi = y[0] + x*x*mu; 
        double p2half = std::abs(energyArg*x*x + y[0] - lambdaArg / (x*x));
        double p      = std::sqrt(2.0*p2half)*Z23;

        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf(phi/(T*x*x));
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        }
        dydx[2] = -std::sqrt(p2half);

        if (x > xmax) {
            double ksi = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])*Z13);
            // std::cout << "x = " << x << ", xmax = " << xmax << ", ksi0 = " << ksi0 << ", ksi = " << ksi << std::endl;
            dydx[3] =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            dydx[3] = -x*x*ksi/p*dydx[3]*dydx[3];
        }

        if (x > xmin && x <= xmax) {

            double ksi2x = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])*Z13); 
            double ksi1x = 1e-8 + std::abs(ksi21 - ksi2x);
            double ax    = ksi2x/ksi21;
            // std::cout << "x = " << x << ", ksi21 = " << ksi21 << ", ksi2x = " << ksi2x << ", ksi1x = " << ksi1x << std::endl;
            
            double Jp13_1 = Jnu(1.0/3.0, ksi1x);
            double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

            double Jp13_2 = Jnu(1.0/3.0, ksi2x);
            double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

            double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
            double R2 = sign_2*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

            if (xmax > 1.0 - 1e-8) {
                dydx[3] = R1;
            }
            else {
                dydx[3]  = (ax*R1 + (1.0 - ax)*R2);
            }

            dydx[3] *= -x*x*dydx[3];
        }

        if (x <= xmin) {
            double ksi = 1e-8 + std::abs(ksi0 + ksi21 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])*Z13);

            // std::cout << "x = " << x << ", ksi = " << ksi << std::endl;

            dydx[3] =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            dydx[3] = -x*x*ksi/p*dydx[3]*dydx[3];
        }

    }
    // rhs at T = 0
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double phi    = y[0] + x*x*mu;
        double p2half = std::abs(energyArg*x*x + y[0] - lambdaArg / (x*x));
        double p      = std::sqrt(2.0*p2half)*Z23;

        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        dydx[2] = -std::sqrt(p2half);

        if (x > xmax) {
            double ksi = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])*Z13);
            dydx[3]    =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            dydx[3]    = -ksi/p*dydx[3]*dydx[3];
        }

        if (x > xmin && x <= xmax) {

            double ksi2x = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])*Z13); 
            double ksi1x = 1e-8 + std::abs(ksi21 - ksi2x);
            double ax    = ksi2x/ksi21;
            
            double Jp13_1 = Jnu(1.0/3.0, ksi1x);
            double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

            double Jp13_2 = Jnu(1.0/3.0, ksi2x);
            double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

            double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
            double R2 = sign_2*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

            if (xmax > 1.0 - 1e-8) {
                dydx[3] = R1;
            }
            else {
                dydx[3]  = (ax*R1 + (1.0 - ax)*R2);
            }

            dydx[3] *= -dydx[3];
        }

        if (x <= xmin) {
            double ksi = 1e-8 + std::abs(ksi0 + ksi21 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])*Z13);
            dydx[3]    =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            dydx[3]    = -ksi/p*dydx[3]*dydx[3];
        }

    }

	double a;
	double energyArg;
	double lambdaArg;
    double mu;
    double V;
    double T;
    double Z;
    double Z13;
    double Z23;
    double r0;
    double ksi0;
    double ksi21;
    double xmin;
    double xmax;
    double sign_2;

    FermiDirac<Half> FDhalf;
    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Knu Knu;
    ::numtk::specfunc::bessel::Knu Ynu;
};

}
}
}
}