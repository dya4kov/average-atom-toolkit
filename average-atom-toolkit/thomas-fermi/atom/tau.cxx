#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/ODE/tau.h>
#include <average-atom-toolkit/thomas-fermi/atom/ODE/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/tau.h>

using namespace aatk::TF;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using ODE::RHSTau;
using ODE::RHSPotential;

double Tau::operator()(const double& _e, const double& _l) {
    setParam(_e, _l);
    if (!ready) setTau();
    return tau;
}

Tau::Tau() :
    V(1.0), T(1.0), Z(1.0), mu(4.100577730112),
    e(1.0), l(1.0),
    tolerance(1e-6) 
{
    ready = false;
}

Tau::Tau(const Tau& t) {
    tolerance = t.tolerance;
    V = t.V; T = t.T;  Z = t.Z; mu = t.mu;
    e = t.e; l = t.l;
    RP = t.RP;
    tau = t.tau; 
    ready  = t.ready;
}

Tau& Tau::operator=(const Tau& t) {
    tolerance = t.tolerance;
    V = t.V; T = t.T; Z = t.Z; mu = t.mu;
    e = t.e; l = t.l;
    RP = t.RP;
    tau = t.tau; 
    ready  = t.ready;
    return *this;
}

void Tau::setV(const double& _V) {
    V = _V;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setV(V);
    ready = false;
}

void Tau::setT(const double& _T) {
    T = _T;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setT(T);
    ready = false;
}

void Tau::setZ(const double& _Z) {
    Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setZ(Z);
    ready = false;
}

void Tau::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    V = _V; T = _T; Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setVTZ(V, T, Z);
    ready = false;
}

void Tau::setTolerance(const double& t) {
    tolerance = t;
    RP.setTolerance(t);
    ready = false;
}

void Tau::setParam(const double& _e, const double& _l) {
    if (std::abs(e - _e) > 1e-10 || std::abs(l - _l)) {
        e = _e; l = _l;
        ready = false; 
    }
}

void Tau::setTau() {
    double xmax = std::sqrt(RP.outer(e, l));
    double xmin = std::sqrt(RP.inner(e, l));
    ready = true;
    tau   = 0.0;

    if (xmax <= xmin) return;

    RHSTau rhsTau;
    Solver<PD853<RHSTau>> tauSolver;
    solver.setTolerance(0.1*tolerance, 0.0);

    RHSPotential rhsPot;
    Solver<PD853<RHSPotential>> potSolver;
    solver.setTolerance(0.1*tolerance, 0.0);

    Array<RHSTau::dim> ytau;
    Array<RHSPotential::dim> ypot;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double r0  = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);

    rhsPot.set_V(V1);
    rhsPot.set_T(T1);
    rhsPot.set_mu(mu1);
    
    rhsTau.set_V(V1);
    rhsTau.set_T(T1);
    rhsTau.set_mu(mu1);
    rhsTau.set_e(e);
    rhsTau.set_l(l);
    
    // split interval
    double alpha = 0.01*(xmax - xmin);

    // (xmin, xmin + alpha)
    ypot.fill(0.0);
    potSolver.setStep(0.1*tolerance);
    potSolver.integrate(rhsPot, ypot, 1.0, xmin);

    double x1    = xmin*xmin;
    double phi1  = ypot[0] + x1*mu1;
    double dphi1 = ypot[1] + mu1;
    double tau1  = r0*std::sqrt(2.0)*alpha/std::sqrt(x1*dphi1 - phi1 + l/x1);

    double tau2  = 0.0;
    double tau3  = 0.0;
    if (xmax < 1.0) {
        // (xmax - alpha, xmax)
        ypot.fill(0.0);
        potSolver.setStep(0.1*tolerance);
        potSolver.integrate(rhsPot, ypot, 1.0, xmax - alpha);
        
        double x2    = xmax*xmax;
        double phi2  = ypot[0] + x2*mu1;
        double dphi2 = ypot[1] + mu1;
        double tau2  = r0*std::sqrt(2.0)*alpha/std::sqrt(x2*dphi2 - phi2 + l/x2);

        // (xmin + alpha, xmax - alpha)
        ytau[0] = ypot[0];
        ytau[1] = ypot[1];
        ytau[2] = 0.0;

        tauSolver.setStep(0.1*tolerance);
        tauSolver.integrate(rhsTau, ytau, xmax - alpha, xmin + alpha);

        tau3 = r0*std::sqrt(2.0)*ytau[2];
    }
    else {
        ytau[0] = 0.0;
        ytau[1] = 0.0;
        ytau[2] = 0.0;

        tauSolver.setStep(0.1*tolerance);
        tauSolver.integrate(rhsTau, ytau, xmax, xmin + alpha);

        tau3 = r0*std::sqrt(2.0)*ytau[2];   
    }

    tau = (tau1 + tau2 + tau3)/Z;
    return;
}