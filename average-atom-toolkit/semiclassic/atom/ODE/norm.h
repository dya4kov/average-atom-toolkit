#pragma once
#include <cmath>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/bessel/Jnu.h>
#include <numeric-toolkit/specfunc/bessel/Ynu.h>
#include <numeric-toolkit/specfunc/bessel/Knu.h>
#include <numeric-toolkit/interpolation/spline.h>

namespace aatk {
namespace semiclassic {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::interpolation::Spline;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSnorm {
	static const Dimension dim    = 4;
    static const Dimension result = 3;

	RHSnorm() : energyArg(1.0), lambdaArg(1.0) {}

    void set_l  (const double& _l)  { lambdaArg = _l; }
    void set_e  (const double& _e)  { energyArg = _e; }
    void set_r0 (const double& _r0) { r0 = _r0; }
    void set_eDens(Spline* _eDens)  { eDens = _eDens; }
    
    void set_ksi0 (const double& _ksi0)  { ksi0  = _ksi0;  }
    void set_ksi21(const double& _ksi21) { ksi21 = _ksi21; }

    void set_RP(const double& _xmin, const double& _xmax) { xmin = _xmin; xmax = _xmax; } 

    void set_sign2(const double& _sign) { sign_2 = _sign; }

	void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { 
        double p2half = std::abs(energyArg*x*x + y[0] - lambdaArg/(x*x));
        double p      = std::sqrt(2.0*p2half);
        auto& density = *eDens;

        dydx[0] = 2.0*x*y[1];
        dydx[1] = x > 0 ? 2.0*density(x)/x : 0.0;
        dydx[2] = -std::sqrt(p2half);

        if (x > xmax) {
            double ksi = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2]));
            // std::cout << "x = " << x << ", xmax = " << xmax << ", ksi0 = " << ksi0 << ", ksi = " << ksi << std::endl;
            dydx[3] =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            dydx[3] = -x*x*ksi/p*dydx[3]*dydx[3];
        }

        if (x > xmin && x <= xmax) {

            double ksi2x = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])); 
            double ksi1x = 1e-8 + std::abs(ksi21 - ksi2x);
            double ax    = ksi2x/ksi21;
            // std::cout << "x = " << x << ", ksi21 = " << ksi21 << ", ksi2x = " << ksi2x << ", ksi1x = " << ksi1x << std::endl;
            
            double Jp13_1 = Jnu(1.0/3.0, ksi1x);
            double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

            double Jp13_2 = Jnu(1.0/3.0, ksi2x);
            double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

            double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
            double R2 = sign_2*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

            dydx[3] = xmax > 1.0 - 1e-8 ? R1 : (ax*R1 + (1.0 - ax)*R2);
            dydx[3] *= -x*x*dydx[3];
        }

        if (x <= xmin) {
            double ksi = 1e-8 + std::abs(ksi0 + ksi21 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2]));

            // std::cout << "x = " << x << ", ksi = " << ksi << std::endl;

            dydx[3] =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            dydx[3] = -x*x*ksi/p*dydx[3]*dydx[3];
        }
    }
private:

	double energyArg;
	double lambdaArg;
    double r0;
    double ksi0;
    double ksi21;
    double xmin;
    double xmax;
    double sign_2;

    Spline* eDens;

    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Knu Knu;
    ::numtk::specfunc::bessel::Knu Ynu;
};

}
}
}