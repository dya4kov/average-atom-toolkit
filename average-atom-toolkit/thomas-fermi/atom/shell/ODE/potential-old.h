#pragma once
#include <cmath>
#include <functional>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>
#include <numeric-toolkit/specfunc/bessel/Jnu.h>
#include <numeric-toolkit/specfunc/bessel/Ynu.h>
#include <numeric-toolkit/specfunc/bessel/Knu.h>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

namespace aatk {
namespace TF {
namespace shell {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::Half;
using ::numtk::specfunc::FD::MHalf;

using ::numtk::specfunc::FermiDiracInc;
using ::numtk::specfunc::FDI::HalfInc2;

using ::aatk::TF::EnergyLevel;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSPotential {
    static const Dimension Nmax = 30;
	static const Dimension dim  = 2 + Nmax*(Nmax + 1)/2 + 2; // 30 levels maximum

	RHSPotential() : V(1.0), T(1.0), Z(1.0), mu(4.100577730112), eb(1.0) {
        rhs = std::bind(&RHSPotential::rhsT, this, _1, _2, _3); eval_a(); eval_r0();
    }

    void set_V   (const double& _V)  { V = _V; eval_a(); eval_r0(); }
    void set_T   (const double& _T)  { T = _T; T12   = std::sqrt(T); T32   = T*T12;
        rhs = T > 1e-10 
                ? std::bind(&RHSPotential::rhsT,  this, _1, _2, _3) 
                : std::bind(&RHSPotential::rhsT0, this, _1, _2, _3);
    }
    void set_Z    (const double& _Z)       { Z = _Z; Z13 = std::pow(Z, 1.0/3.0); Z23 = Z13*Z13; Z43 = Z23*Z23; }
    void set_mu   (const double& _mu)      { mu = _mu; }
    void set_eb   (const double& _eb)      { eb = _eb; }
    void set_nmax (const int& _nmax) { nmax = Nmax > _nmax ? _nmax : Nmax; }
    
    void set_e    (const std::vector<double>& _e)     { e = _e; }
    void set_ksi0 (const std::vector<double>& _ksi0)  { ksi0  = _ksi0;  }
    void set_ksi21(const std::vector<double>& _ksi21) { ksi21 = _ksi21; }

    void set_RP   (const std::vector<double>& _rpi, 
                   const std::vector<double>& _rpo) { 
        rpi = _rpi; 
        rpo = _rpo; 
    } 

    void set_sign (const std::vector<double>& _sign)  { sign = _sign; }

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
        
        double phi   = y[0] + x2*mu; 
        double phiBE = y[0] + x2*eb;

        dydx.fill(0.0);

        dydx[0] = 2.0*x*y[1];
        if (x > 0) {
            dydx[1] = 2.0*a*x3*T32*FDhalf(phi/(T*x2));
        }
        else {
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
        }

        double rhoSC = 0.0;
        for (int n = 1; n < nmax; ++n) {
            for (int l = 0; l < n; ++l) {

                int i = 2 + n*(n - 1)/2 + l;
                if (e[i] > eb) continue;

                double lambda = 0.5*(l + 0.5)*(l + 0.5) / r0 / r0 / Z23;
                double p2half = std::abs(e[i]*x2 + y[0] - lambda / (x2));
                double p      = std::sqrt(2.0*p2half)*Z23;
                dydx[2 + i]   = -std::sqrt(p2half);
                double Nnl    = 2.0*(2.0*l + 1.0)/(1.0 + std::exp((e[i] - mu)/T));
                double Rnl    = R(n, l, p, x, y);
                rhoSC        += Nnl*Rnl*Rnl/x;
            }
        }

        double rhoTF = 0.0;
        if (phiBE <= 0.0) rhoTF = 0.0;
        else {
            double argX = phi/(T*x2);
            double argY = phiBE/(T*x2);
            if (argY - argX > 25.0)
                rhoTF = a*x3*T32*FDhalf(argX)*Z*Z;
            else
                rhoTF = a*x3*T32*FDhalfI(argX, argY)*Z*Z;
        }

        std::cout << "x = " << x <<  ", rhoTF = " << rhoTF << ", rhoSC = " << rhoSC << std::endl;

        double drho_sh = rhoSC - rhoTF;
        double drho_tf = a*x*T12*FDmhalf(phi/(T*x2))*Z23;

        dydx[dim - 2] = 2.0*x*y[dim - 1];
        dydx[dim - 1] = drho_tf*y[dim - 2] + 2.0*drho_sh;

        dydx[dim - 2] = 0.0;
        dydx[dim - 1] = 0.0;

        // std::cout << "x = " << x << ", y0 = " << y[0] << ", y1 = " << y[1] << ", y[2] = " << y[2] << ", yn = " << y[dim - 1] << std::endl;
    }
    // rhs at T = 0
    void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
        double x2    = x*x;
        double x3    = x*x2;
        double x5    = x2*x3;
        
        double phi   = y[0] + x2*mu; 
        double phiBE = y[0] + x2*eb;

        dydx.fill(0.0);

        dydx[0] = 2.0*x*y[1];
        dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;

        double rhoSC = 0.0;
        for (int n = 1; n < nmax; ++n) {
            for (int l = 0; l < n; ++l) {

                int i = 2 + n*(n - 1)/2 + l;
                if (e[i] > eb) continue;

                double lambda = 0.5*(l + 0.5)*(l + 0.5) / r0 / r0 / Z23;
                double p2half = std::abs(e[i]*x2 + y[0] - lambda / (x2));
                double p      = std::sqrt(2.0*p2half)*Z23;
                dydx[2 + i]   = -std::sqrt(p2half);
                double Nnl    = 2.0*(2.0*l + 1.0);
                double Rnl    = R(n, l, p, x, y);
                rhoSC        += Nnl*Rnl*Rnl/x;
            }
        }

        double rhoTF = 0.0;

        if (phiBE <= 0.0) rhoTF = 0.0;
        else {
            rhoTF = 2.0/3.0*a*x3*Z*Z*std::pow(std::min(eb, mu)*x2 + y[0], 1.5);
        }

        double drho_sh = rhoSC - rhoTF;
        double drho_tf = 2.0*a*std::sqrt(phi)*Z23;

        dydx[dim - 2] = 2.0*x*y[dim - 1];
        dydx[dim - 1] = drho_tf*y[dim - 2] + 2.0*drho_sh;

    }

    double R(int n, int l, double p, double x, Array<dim>& y) {

        int i = 2 + n*(n - 1)/2 + l;
        double result = 0.0;

        if (x > rpo[i]) {
            double ksi = 1e-8 + std::abs(ksi0[i] - 2.0*r0*std::sqrt(2.0)*std::abs(y[i])*Z13);
            result     = ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            result     = ksi/p*result;
        }

        if (x > rpi[i] && x <= rpo[i]) {

            double ksi2x = 1e-8 + std::abs(ksi0[i] - 2.0*r0*std::sqrt(2.0)*std::abs(y[i])*Z13); 
            double ksi1x = 1e-8 + std::abs(ksi21[i] - ksi2x);
            double ax    = ksi2x/ksi21[i];
            
            double Jp13_1 = Jnu(1.0/3.0, ksi1x);
            double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

            double Jp13_2 = Jnu(1.0/3.0, ksi2x);
            double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

            double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
            double R2 = sign[i]*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

            result = (ax*R1 + (1.0 - ax)*R2);
        }

        if (x <= rpi[i]) {
            double ksi = 1e-8 + std::abs(ksi0[i] + ksi21[i] - 2.0*r0*std::sqrt(2.0)*std::abs(y[i])*Z13);
            result     = ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
            result     = ksi/p*result;
        }

        return result;
    }

    int nmax;

	double a;
    double mu;
    double eb;
    double V;
    double T, T12, T32;
    double Z, Z13, Z23, Z43;
    double r0;

    std::vector<double>  ksi0;
    std::vector<double> ksi21;
    std::vector<double>  sign;
    std::vector<double>   rpi;
    std::vector<double>   rpo;
    std::vector<double>     e;

    FermiDirac<Half>  FDhalf;
    FermiDirac<MHalf> FDmhalf;
    FermiDiracInc<HalfInc2> FDhalfI;
    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Knu Knu;
    ::numtk::specfunc::bessel::Knu Ynu;

};

}
}
}
}