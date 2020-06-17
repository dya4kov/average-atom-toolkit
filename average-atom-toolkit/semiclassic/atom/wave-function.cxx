#include <cmath>
#include <numeric>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/semiclassic/atom.h>
#include <average-atom-toolkit/semiclassic/atom/ODE/ksi.h>
#include <average-atom-toolkit/semiclassic/atom/ODE/norm.h>

namespace aatk {
namespace semiclassic {

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using aatk::semiclassic::ODE::RHSksi;
using aatk::semiclassic::ODE::RHSnorm;

std::vector<double> Atom::waveFunction(double energy, double lambda, const std::vector<double>& x) {
    std::vector<double> result(x.size(), 0.0);
	waveFunction(energy, lambda, x.data(), result.data(), x.size());
    return result;
}
std::vector<double> Atom::waveFunctionVec(double energy, double lambda, const std::vector<double>& x) {
    std::vector<double> result(x.size(), 0.0);
    waveFunction(energy, lambda, x.data(), result.data(), x.size());
    return result;
}
double Atom::waveFunction(double energy, double lambda, double x) {
    double result;
    waveFunction(energy, lambda, &x, &result, 1);
	return result;
}
void Atom::waveFunction(double energy, double lambda, const double* xin, double* result, std::size_t n) {

    auto rpo = outerRP(energy, lambda);
    auto rpi = innerRP(energy, lambda);

    double normValue = waveFunctionNorm(
    	energy, lambda, 
    	rpi[0], rpo[0],
    	rpo[1], rpo[2]
    );

    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double> x(n);

    std::sort(idx.begin(), idx.end(),
       [xin](std::size_t i1, std::size_t i2) {
        return xin[i1] > xin[i2]; 
       }
    );

    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::sqrt(xin[idx[i]]);
        result[i] = 0.0;
    }

    double xmax = std::sqrt(rpo[0]);
    double xmin = std::sqrt(rpi[0]);

    if (xmax <= xmin) return;

    double l = 0.5*lambda*lambda / r0 / r0;
    double e = energy;

    RHSksi rhs;
    rhs.set_e(e);
    rhs.set_l(l);
    rhs.set_eDens(densityInterpolation);

    Solver<PD853<RHSksi>> solver;
    solver.setTolerance(0.1*tolerance, 0.0);

    Array<RHSksi::dim> yksi;

    double ksi0 = 0.0;

    if (xmax < 1.0 - tolerance) {
        yksi[0] = 0.0; 
        yksi[1] = 0.0;
        yksi[2] = 0.0;
        
        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, 1.0, xmax);

        ksi0 = 2.0*r0*std::sqrt(2.0)*yksi[2];
    }

    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Ynu Ynu;
    ::numtk::specfunc::bessel::Knu Knu;

    yksi[0] = rpo[1];
    yksi[1] = rpo[2];
    yksi[2] = 0.0;

    solver.setStep(tolerance);
    solver.integrate(rhs, yksi, xmax, xmin);

    double ksi21 = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);
    double Jp13 = Jnu(1.0/3.0, ksi21);
    double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21));
    double sign_2 = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;

    // gamma(1/3)
    double gamma   = 2.6789385347077476336557;

    double RPOx    = xmax*xmax;
    double RPOphi  = rpo[1];
    double RPOdphi = rpo[2];
    double RPOc    = 2.0/RPOx*((RPOdphi - RPOphi/RPOx) + 2.0*l/(RPOx*RPOx));

    double ksi16pm12_2 = std::pow(2.0*r0/(3.0*std::abs(RPOc)), 1.0/6.0);

    double R2approx = sign_2*normValue*gamma/M_PI*ksi16pm12_2*std::pow(2.0, -2.0/3.0);

    double RPIx    = xmin*xmin;
    double RPIphi  = rpi[1];
    double RPIdphi = rpi[2];
    double RPIc    = 2.0/RPIx*((RPIdphi - RPIphi/RPIx) + 2.0*l/(RPIx*RPIx));

    double ksi16pm12_1 = std::pow(2.0*r0/(3.0*std::abs(RPIc)), 1.0/6.0);

    double R1approx = normValue*gamma/M_PI*ksi16pm12_1*std::pow(2.0, -2.0/3.0);

    double ksiParam = std::pow(2.0, 1.0/3.0)*std::sqrt(3.0)*M_PI/(gamma*gamma);
    double ksiRPeps = 2e-3;

    double xFrom = 1.0;
    std::size_t i = 0;

    yksi[0] = 0.0; 
    yksi[1] = 0.0;
    yksi[2] = 0.0;

    while (x[i] > xmax && i < n) {
        double xTo = x[i];

        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, xFrom, xTo);

        double ksi = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*yksi[2]);
        double p = std::sqrt(2.0*std::abs(e*xTo*xTo + yksi[0] - l / (xTo*xTo)))/xTo;

        if      (ksi > 100.0)    result[idx[i]] = 0.0;
        else if (ksi < ksiRPeps) result[idx[i]] = R2approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
        else    result[idx[i]] = sign_2*normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);

        ++i; xFrom = xTo;
    }

    yksi[0] = rpo[1];
    yksi[1] = rpo[2];
    yksi[2] = 0.0;

    xFrom = xmax;

    while (x[i] > xmin && x[i] <= xmax && i < n) {

        double xTo = x[i];

        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, xFrom, xTo);

        double ksi2x  = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]); 
        double a = ksi2x/ksi21;

        double p = std::sqrt(2.0*std::abs(e*xTo*xTo + yksi[0] - l / (xTo*xTo)))/xTo;

        double ksi1x  = 1e-8 + std::abs(ksi21 - ksi2x);

        double Jp13_1 = Jnu(1.0/3.0, ksi1x);
        double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

        double Jp13_2 = Jnu(1.0/3.0, ksi2x);
        double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

        double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
        double R2 = sign_2*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

        result[idx[i]] = normValue*(a*R1 + (1.0 - a)*R2);

        if (xmax > 1.0 - tolerance) {
            result[idx[i]] = normValue*R1;
        }

        if (ksi2x < ksiRPeps && xmax < 1.0 - tolerance) 
            result[idx[i]] = R2approx*(1.0 + ksiParam*std::pow(ksi2x, 2.0/3.0));

        if (ksi1x < ksiRPeps) 
            result[idx[i]] = R1approx*(1.0 + ksiParam*std::pow(ksi1x, 2.0/3.0));

        ++i; xFrom = xTo;
    }

    yksi[0] = rpi[1];
    yksi[1] = rpi[2];
    yksi[2] = 0.0;

    xFrom = xmin;

    while (x[i] < xmin && i < n) {

        double xTo = x[i];

        if (x[i]*x[i] > tolerance) {

            solver.setStep(tolerance);
            solver.integrate(rhs, yksi, xFrom, xTo);

            double ksi = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);
            double p = std::sqrt(2.0*std::abs(e*xTo*xTo + yksi[0] - l / (xTo*xTo)))/xTo;

            if      (ksi > 100.0)    result[idx[i]] = 0.0;
            else if (ksi < ksiRPeps) result[idx[i]] = R1approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
            else    result[idx[i]] = normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);
        }
        else { result[idx[i]] = 0.0; }

        ++i; xFrom = xTo;
    }
}

double Atom::waveFunctionNorm(
	double energy, double lambda, 
	double rpi,    double rpo, 
	double potrpo, double dpotrpo
) {
	double xmax = std::sqrt(rpo);
    double xmin = std::sqrt(rpi);

    if (xmax <= xmin) return 1.0;

    double l = 0.5*lambda*lambda / r0 / r0;
    double e = energy;

    RHSksi rhs;
    rhs.set_e(e);
    rhs.set_l(l);
    rhs.set_eDens(densityInterpolation);

    Solver<PD853<RHSksi>> solver;
    solver.setTolerance(0.1*tolerance, 0.0);

    Array<RHSksi::dim> yksi;

    double ksi0 = 0.0;

    if (xmax < 1.0 - tolerance) {
        yksi[0] = 0.0; 
        yksi[1] = 0.0;
        yksi[2] = 0.0;
        
        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, 1.0, xmax);

        ksi0 = 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);
    }


    yksi[0] = potrpo;
    yksi[1] = dpotrpo;
    yksi[2] = 0.0;

    solver.setStep(tolerance);
    solver.integrate(rhs, yksi, xmax, xmin);


    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Ynu Ynu;
    ::numtk::specfunc::bessel::Knu Knu;

    double ksi21 = 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);
    double Jp13 = Jnu(1.0/3.0, ksi21);
    double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21));
    double sign_2 = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;


    RHSnorm rhsNorm;
    rhsNorm.set_e(e);
    rhsNorm.set_l(l);
    rhsNorm.set_r0(r0);
    rhsNorm.set_eDens(densityInterpolation);
    rhsNorm.set_ksi0(ksi0);
    rhsNorm.set_ksi21(ksi21);
    rhsNorm.set_RP(xmin, xmax);
    rhsNorm.set_sign2(sign_2);

    Solver<PD853<RHSnorm>> normSolver;
    normSolver.setTolerance(0.1*tolerance, 0.0);


    Array<RHSnorm::dim> ynorm;
    ynorm.fill(0.0);

    normSolver.setStep(tolerance);
    normSolver.integrate(rhsNorm, ynorm, 1.0, std::sqrt(0.1*tolerance));

    return std::sqrt(1.0/(4.0*r0*std::abs(ynorm[3])));
}

}
}