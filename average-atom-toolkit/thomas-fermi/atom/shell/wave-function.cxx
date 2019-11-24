#include <cmath>
#include <numeric>
#include <algorithm>

#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/bessel/Jnu.h>
#include <numeric-toolkit/specfunc/bessel/Knu.h>
#include <numeric-toolkit/specfunc/bessel/Ynu.h>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/shell/wave-function.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/ODE/ksi.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/ODE/norm.h>

using namespace aatk::TF::shell;

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using numtk::specfunc::FermiDirac;
using numtk::specfunc::FD::Half;
using aatk::TF::shell::ODE::RHSksi;
using aatk::TF::shell::ODE::RHSnorm;

WaveFunction::WaveFunction() :
    V(1.0), T(1.0), Z(1.0), mu(4.100577730112),
    energy(1.0), lambda(1.0),
    tolerance(1e-6) 
{
    ready = false;
}

WaveFunction::WaveFunction(const WaveFunction& wf) {
    tolerance = wf.tolerance;
    V = wf.V; T = wf.T;  Z = wf.Z; mu = wf.mu;
    energy = wf.energy; lambda = wf.lambda;
    RP = wf.RP;
    normValue = wf.normValue; 
    ready = wf.ready;
}

WaveFunction& WaveFunction::operator=(const WaveFunction& wf) {
    tolerance = wf.tolerance;
    V = wf.V; T = wf.T; Z = wf.Z; mu = wf.mu;
    energy = wf.energy; lambda = wf.lambda;
    RP = wf.RP;
    normValue = wf.normValue; 
    ready = wf.ready;
    return *this;
}

void WaveFunction::setV(const double _V) {
    V = _V;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setV(V);
    ready = false;
}

void WaveFunction::setT(const double _T) {
    T = _T;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setT(T);
    ready = false;
}

void WaveFunction::setZ(const double _Z) {
    Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setZ(Z);
    ready = false;
}

void WaveFunction::setVTZ(
    const double _V,
    const double _T,
    const double _Z
) {
    V = _V; T = _T; Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setVTZ(V, T, Z);
    ready = false;
}

void WaveFunction::setTolerance(const double t) {
    tolerance = t;
    RP.setTolerance(t);
    ready = false;
}

void WaveFunction::setParam(const double _energy, const double _lambda) {
    if (std::abs(energy - _energy) > 1e-10 || std::abs(lambda - _lambda) > 1e-10) {
        energy = _energy; lambda = _lambda;
        ready = false; 
    }
}

double WaveFunction::norm(const double _energy, const double _lambda) {
    setParam(_energy, _lambda);
    if (!ready) setNorm();
    return normValue;
}

void WaveFunction::setNorm() {

    double xmax = std::sqrt(RP.outer(energy, lambda));
    double xmin = std::sqrt(RP.inner(energy, lambda));

    normValue = 1.0;
    if (xmax <= xmin) return;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);
    double  r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double   l = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double   e = energy*std::pow(Z, -4.0/3.0);

    RHSksi rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    rhs.set_e(e);
    rhs.set_l(l);

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

        ksi0 = 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0);
    }

    yksi[0] = RP.outerY(energy, lambda)[0];
    yksi[1] = RP.outerY(energy, lambda)[1];
    yksi[2] = 0.0;

    solver.setStep(tolerance);
    solver.integrate(rhs, yksi, xmax, xmin);

    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Ynu Ynu;
    ::numtk::specfunc::bessel::Knu Knu;

    double ksi21 = 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0);
    double Jp13 = Jnu(1.0/3.0, ksi21);
    double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21));
    double sign_2 = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;

    // std::cout << "ksi0 = " << ksi0 << ", ksi21 = " << ksi21 << ", xmin = " << xmin << ", xmax = " << xmax << std::endl;

    RHSnorm rhsNorm;
    rhsNorm.set_V(V1);
    rhsNorm.set_T(T1);
    rhsNorm.set_Z(Z);
    rhsNorm.set_mu(mu1);
    rhsNorm.set_e(e);
    rhsNorm.set_l(l);
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

    r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
    normValue = std::sqrt(1.0/(4.0*r0*std::abs(ynorm[3])));

    ready = true;

    // std::cout << "norm = " << normValue << std::endl; 
}

double WaveFunction::operator()(const double _energy, const double _lambda, const double _x) {
    if (_x < tolerance) return 0.0;
    setParam(_energy, _lambda);
    if (!ready) setNorm();

    double x    = std::sqrt(_x);
    double xmax = std::sqrt(RP.outer(energy, lambda));
    double xmin = std::sqrt(RP.inner(energy, lambda));

    if (xmax <= xmin) return 0.0;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);
    double  r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double   l = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double   e = energy*std::pow(Z, -4.0/3.0);

    double ksi0 = 0.0;
    double result = 0.0;

    RHSksi rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    rhs.set_e(e);
    rhs.set_l(l);

    Solver<PD853<RHSksi>> solver;
    solver.setTolerance(0.1*tolerance, 0.0);

    Array<RHSksi::dim> yksi;

    if (xmax < 1.0 - tolerance) {
        yksi[0] = 0.0; 
        yksi[1] = 0.0;
        yksi[2] = 0.0;
        
        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, 1.0, xmax);

        ksi0 = 2.0*r0*std::sqrt(2.0)*yksi[2]*std::pow(Z, 1.0/3.0);
    }

    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Ynu Ynu;
    ::numtk::specfunc::bessel::Knu Knu;

    yksi[0] = RP.outerY(energy, lambda)[0];
    yksi[1] = RP.outerY(energy, lambda)[1];
    yksi[2] = 0.0;

    solver.setStep(tolerance);
    solver.integrate(rhs, yksi, xmax, xmin);

    double ksi21 = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0);
    double Jp13 = Jnu(1.0/3.0, ksi21);
    double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21));
    double sign_2 = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;

    // gamma(1/3)
    double gamma   = 2.6789385347077476336557;

    double RPOx    = xmax*xmax;
    double RPOphi  = RP.outerY(energy, lambda)[0] + xmax*xmax*mu1;
    double RPOdphi = RP.outerY(energy, lambda)[1] + mu1;
    double RPOc    = 2.0/RPOx*((RPOdphi - RPOphi/RPOx) + 2.0*l/(RPOx*RPOx));

    double ksi16pm12_2 = std::pow(2.0*r0/(3.0*std::abs(RPOc)), 1.0/6.0)*std::pow(Z, -5.0/18.0);

    double R2approx = sign_2*normValue*gamma/M_PI*ksi16pm12_2*std::pow(2.0, -2.0/3.0);

    double RPIx    = xmin*xmin;
    double RPIphi  = RP.innerY(energy, lambda)[0] + xmin*xmin*mu1;
    double RPIdphi = RP.innerY(energy, lambda)[1] + mu1;
    double RPIc    = 2.0/RPIx*((RPIdphi - RPIphi/RPIx) + 2.0*l/(RPIx*RPIx));

    double ksi16pm12_1 = std::pow(2.0*r0/(3.0*std::abs(RPIc)), 1.0/6.0)*std::pow(Z, -5.0/18.0);;

    double R1approx = normValue*gamma/M_PI*ksi16pm12_1*std::pow(2.0, -2.0/3.0);

    double ksiParam = std::pow(2.0, 1.0/3.0)*std::sqrt(3.0)*M_PI/(gamma*gamma);
    double ksiRPeps = 2e-3;

    if (x > xmax) {
        yksi[0] = 0.0; 
        yksi[1] = 0.0;
        yksi[2] = 0.0;

        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, 1.0, x);

        double ksi = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0));
        double p = std::sqrt(2.0*std::abs(e*x*x + yksi[0] - l / (x*x)))/x*std::pow(Z, 2.0/3.0);

        if      (ksi > 100.0)    result = 0.0;
        else if (ksi < ksiRPeps) result = R2approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
        else    result = sign_2*normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);

    }

    if (x > xmin && x <= xmax) {

        yksi[0] = RP.outerY(energy, lambda)[0];
        yksi[1] = RP.outerY(energy, lambda)[1];
        yksi[2] = 0.0;

        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, xmax, x);

        double ksi2x = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0); 
        double ksi1x = 1e-8 + std::abs(ksi21 - ksi2x);

        double a = ksi2x/ksi21;
        double p = std::sqrt(2.0*std::abs(e*x*x + yksi[0] - l / (x*x)))/x*std::pow(Z, 2.0/3.0);

        double Jp13_1 = Jnu(1.0/3.0, ksi1x);
        double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

        double Jp13_2 = Jnu(1.0/3.0, ksi2x);
        double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

        double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
        double R2 = sign_2*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

        result = normValue*(a*R1 + (1.0 - a)*R2);

        if (xmax > 1.0 - tolerance) {
            result = normValue*R1;
        }

        if (ksi2x < ksiRPeps && xmax < 1.0 - tolerance) 
            result = R2approx*(1.0 + ksiParam*std::pow(ksi2x, 2.0/3.0));

        if (ksi1x < ksiRPeps) 
            result = R1approx*(1.0 + ksiParam*std::pow(ksi1x, 2.0/3.0));
    }

    if (x < xmin) {
        yksi[0] = RP.innerY(energy, lambda)[0];
        yksi[1] = RP.innerY(energy, lambda)[1];
        yksi[2] = 0.0;

        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, xmin, x);

        double ksi = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0);
        double p = std::sqrt(2.0*std::abs(e*x*x + yksi[0] - l / (x*x)))/x*std::pow(Z, 2.0/3.0);

        if      (ksi > 100.0)    result = 0.0;
        else if (ksi < ksiRPeps) result = R1approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
        else    result = normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);
    }

    return result;
}

void WaveFunction::operator()(const double _energy, const double _lambda, const double* _x, double* result, const std::size_t n) {
    setParam(_energy, _lambda);
    if (!ready) setNorm();
    
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* x      = new double[n];

    std::sort(idx.begin(), idx.end(),
       [_x](std::size_t i1, std::size_t i2) {
        return _x[i1] > _x[i2]; 
       }
    );

    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::sqrt(_x[idx[i]]);
        result[i] = 0.0;
    }

    double xmax = std::sqrt(RP.outer(energy, lambda));
    double xmin = std::sqrt(RP.inner(energy, lambda));

    if (xmax <= xmin) return;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);
    double  r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double   l = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double   e = energy*std::pow(Z, -4.0/3.0);

    RHSksi rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    rhs.set_e(e);
    rhs.set_l(l);

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

        ksi0 = 2.0*r0*std::sqrt(2.0)*yksi[2]*std::pow(Z, 1.0/3.0);
    }

    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Ynu Ynu;
    ::numtk::specfunc::bessel::Knu Knu;

    yksi[0] = RP.outerY(energy, lambda)[0];
    yksi[1] = RP.outerY(energy, lambda)[1];
    yksi[2] = 0.0;

    solver.setStep(tolerance);
    solver.integrate(rhs, yksi, xmax, xmin);

    double ksi21 = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0);
    double Jp13 = Jnu(1.0/3.0, ksi21);
    double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21));
    double sign_2 = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;

    // gamma(1/3)
    double gamma   = 2.6789385347077476336557;

    double RPOx    = xmax*xmax;
    double RPOphi  = RP.outerY(energy, lambda)[0] + xmax*xmax*mu1;
    double RPOdphi = RP.outerY(energy, lambda)[1] + mu1;
    double RPOc    = 2.0/RPOx*((RPOdphi - RPOphi/RPOx) + 2.0*l/(RPOx*RPOx));

    double ksi16pm12_2 = std::pow(2.0*r0/(3.0*std::abs(RPOc)), 1.0/6.0)*std::pow(Z, -5.0/18.0);

    double R2approx = sign_2*normValue*gamma/M_PI*ksi16pm12_2*std::pow(2.0, -2.0/3.0);

    double RPIx    = xmin*xmin;
    double RPIphi  = RP.innerY(energy, lambda)[0] + xmin*xmin*mu1;
    double RPIdphi = RP.innerY(energy, lambda)[1] + mu1;
    double RPIc    = 2.0/RPIx*((RPIdphi - RPIphi/RPIx) + 2.0*l/(RPIx*RPIx));

    double ksi16pm12_1 = std::pow(2.0*r0/(3.0*std::abs(RPIc)), 1.0/6.0)*std::pow(Z, -5.0/18.0);;

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

        double ksi = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*yksi[2]*std::pow(Z, 1.0/3.0));
        double p = std::sqrt(2.0*std::abs(e*xTo*xTo + yksi[0] - l / (xTo*xTo)))/xTo*std::pow(Z, 2.0/3.0);

        if      (ksi > 100.0)    result[idx[i]] = 0.0;
        else if (ksi < ksiRPeps) result[idx[i]] = R2approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
        else    result[idx[i]] = sign_2*normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);

        ++i; xFrom = xTo;
    }

    yksi[0] = RP.outerY(energy, lambda)[0];
    yksi[1] = RP.outerY(energy, lambda)[1];
    yksi[2] = 0.0;

    xFrom = xmax;

    while (x[i] > xmin && x[i] <= xmax && i < n) {

        double xTo = x[i];

        solver.setStep(tolerance);
        solver.integrate(rhs, yksi, xFrom, xTo);

        double ksi2x  = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0); 
        double a = ksi2x/ksi21;

        double p = std::sqrt(2.0*std::abs(e*xTo*xTo + yksi[0] - l / (xTo*xTo)))/xTo*std::pow(Z, 2.0/3.0);

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

    yksi[0] = RP.innerY(energy, lambda)[0];
    yksi[1] = RP.innerY(energy, lambda)[1];
    yksi[2] = 0.0;

    xFrom = xmin;

    while (x[i] < xmin && i < n) {

        double xTo = x[i];

        if (x[i]*x[i] > tolerance) {

            solver.setStep(tolerance);
            solver.integrate(rhs, yksi, xFrom, xTo);

            double ksi = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0);
            double p = std::sqrt(2.0*std::abs(e*xTo*xTo + yksi[0] - l / (xTo*xTo)))/xTo*std::pow(Z, 2.0/3.0);

            if      (ksi > 100.0)    result[idx[i]] = 0.0;
            else if (ksi < ksiRPeps) result[idx[i]] = R1approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
            else    result[idx[i]] = normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);
        }
        else { result[idx[i]] = 0.0; }

        ++i; xFrom = xTo;
    }
}

double* WaveFunction::operator()(const double _energy, const double _lambda, const double* _x, const std::size_t n) {
    double* result = new double[n];
    operator()(_energy, _lambda, _x, result, n);
    return result;
}

std::vector<double> WaveFunction::operator()(const double _energy, const double _lambda, const std::vector<double>& _x) {
    std::vector<double> result(_x.size());
    operator()(_energy, _lambda, _x.data(), result.data(), _x.size());
    return result;
}
