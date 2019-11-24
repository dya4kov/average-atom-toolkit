#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/shell/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/ODE/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/electron-density.h>
#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/shell/chemical-potential.h>

using namespace aatk::TF::shell;
using ::numtk::ODE::Array;
using ::numtk::ODE::Solver;
using ::numtk::ODE::stepper::PD853;

using ::aatk::TF::shell::ODE::RHSPotential;

#define DRHOSH_INTERPOLATION_SIZE 3501

Potential::Potential() : 
    V(1.0), T(1.0), Z(1.0), nmax(25), 
    tolerance(1e-6), ready(false), 
    x(DRHOSH_INTERPOLATION_SIZE, 0.0),
    sqrtx(DRHOSH_INTERPOLATION_SIZE, 0.0), 
    drhosh(DRHOSH_INTERPOLATION_SIZE, 0.0)
{
	double dx = 1.0/(x.size() - 1.0);
	for (std::size_t i = 0; i < x.size(); ++i) {
		sqrtx[i] = i*dx;
		x[i] = sqrtx[i]*sqrtx[i];
	}
}

Potential::Potential(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    mu = p.mu; 
    dmush = p.dmush;
    nmax = p.nmax;
    tolerance = p.tolerance;
    ready = p.ready;
}

Potential& Potential::operator=(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    mu = p.mu; 
    dmush = p.dmush;
    nmax = p.nmax;
    tolerance = p.tolerance;
    ready = p.ready;
    return *this;
}

void Potential::setV(const double& _V) { 
    if (std::abs(V - _V) > tolerance) {
        V = _V; ready = false;
    }
}

void Potential::setT(const double& _T) {
    if (std::abs(T - _T) > tolerance) {
        T = _T; ready = false;
    }
}

void Potential::setZ(const double& _Z) {
    if (std::abs(Z - _Z) > tolerance) {
        Z = _Z; ready = false;
    }
}

void Potential::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    if (std::abs(V - _V) > tolerance ||
        std::abs(T - _T) > tolerance ||
        std::abs(Z - _Z) > tolerance) {
        V = _V; T = _T; Z = _Z;
        ready = false;
    }
}

void Potential::setTolerance(const double& t) {
    if (std::abs(tolerance - t) > 1e-15) {
        tolerance = std::max(t, tolerance);
        ready = false;
    }
}

void Potential::setNmax(const int& _nmax) {
    if (_nmax != nmax) {
        nmax = _nmax; 
        ready = false;
    }
}

void Potential::prepare() {

    // prepare electron states
    ::aatk::TF::ElectronStates N;
    N.setVTZ(V, T, Z);
    N.setTolerance(tolerance);
    N.setNmax(nmax);

    eb = N.eBoundary();
    // chemical potential

    ::aatk::TF::ChemicalPotential M;
    M.setTolerance(tolerance);
    M.setZ(Z);
    mu = M(V, T);

    ::aatk::TF::shell::ChemicalPotential dM;
    dM.setTolerance(tolerance);
    dM.setZ(Z);
    dM.setNmax(nmax);
    dmush = dM(V, T);

    // prepare electron density
    auto eLevel = N.eLevel();
    ::aatk::TF::shell::ElectronDensity drho;
    drho.setVTZ(V, T, Z);
    drho.setTolerance(tolerance);
    drho.setNmax(nmax);
    drho.setEnergyLevels(N.eLevel());
    drho.setBoundary(eb);

    double r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);

    drhosh = drho(x);
    drhosh[0] = 0.0;
    for (std::size_t i = 1; i < drhosh.size(); ++i) {
    	drhosh[i] = 4.0*M_PI*r0*r0*x[i]*x[i]*drhosh[i];
    }
    ready = true;
}

double Potential::operator()(const double& x) {
    if (!ready) prepare();

    RHSPotential rhs;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V    (V1);
    rhs.set_T    (T1);
    rhs.set_Z    (Z);
    rhs.set_mu   (mu1);
    rhs.set_eb   (eb*std::pow(Z, -4.0/3.0));
    rhs.set_drho (sqrtx, drhosh);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);

    psi.fill(0.0);
    psi[2] = psi[3] = dmush;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, tolerance);

    if (xFrom > xTo) solver.integrate(rhs, psi, xFrom, xTo);
    return psi[2]/x - dmush;
}

double Potential::dx(const double& x) {
    if (!ready) prepare();

    RHSPotential rhs;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V    (V1);
    rhs.set_T    (T1);
    rhs.set_Z    (Z);
    rhs.set_mu   (mu1);
    rhs.set_eb   (eb*std::pow(Z, -4.0/3.0));
    rhs.set_drho (sqrtx, drhosh);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);

    psi.fill(0.0);
    psi[2] = psi[3] = dmush;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, tolerance);

    if (xFrom > xTo) solver.integrate(rhs, psi, xFrom, xTo);
    return psi[3]/x - dmush;
}

double* Potential::operator()(const double* x, const std::size_t& n) {

    if (!ready) prepare();

    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    RHSPotential rhs;

    rhs.set_V    (V1);
    rhs.set_T    (T1);
    rhs.set_Z    (Z);
    rhs.set_mu   (mu1);
    rhs.set_eb   (eb*std::pow(Z, -4.0/3.0));
    rhs.set_drho (sqrtx, drhosh);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;

    psi.fill(0.0);
    psi[2] = psi[3] = dmush;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        result[i] = psi[2]/x[i] - dmush;
        xFrom = xTo;
    }

    return result;
}

double* Potential::dx(const double* x, const std::size_t& n) {

    if (!ready) prepare();

    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    RHSPotential rhs;

    rhs.set_V    (V1);
    rhs.set_T    (T1);
    rhs.set_Z    (Z);
    rhs.set_mu   (mu1);
    rhs.set_eb   (eb*std::pow(Z, -4.0/3.0));
    rhs.set_drho (sqrtx, drhosh);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;

    psi.fill(0.0);
    psi[2] = psi[3] = dmush;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        result[i] = psi[3]/x[i] - dmush;
        xFrom = xTo;
    }

    return result;
}

std::vector<double> Potential::operator()(const std::vector<double>& x) {

    if (!ready) prepare();

    std::size_t n = x.size();
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double> result(n);

    std::sort(idx.begin(), idx.end(),
       [&x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    RHSPotential rhs;

    rhs.set_V    (V1);
    rhs.set_T    (T1);
    rhs.set_Z    (Z);
    rhs.set_mu   (mu1);
    rhs.set_eb   (eb*std::pow(Z, -4.0/3.0));
    rhs.set_drho (sqrtx, drhosh);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;

    psi.fill(0.0);
    psi[2] = psi[3] = dmush;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        result[i] = psi[2]/x[i] - dmush;
        xFrom = xTo;
    }

    return result;
}

std::vector<double> Potential::dx(const std::vector<double>& x) {

    if (!ready) prepare();

    std::size_t n = x.size();
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double> result(n);

    std::sort(idx.begin(), idx.end(),
       [&x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    RHSPotential rhs;

    rhs.set_V    (V1);
    rhs.set_T    (T1);
    rhs.set_Z    (Z);
    rhs.set_mu   (mu1);
    rhs.set_eb   (eb*std::pow(Z, -4.0/3.0));
    rhs.set_drho (sqrtx, drhosh);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;

    psi.fill(0.0);
    psi[2] = psi[3] = dmush;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        result[i] = psi[3]/x[i] - dmush;
        xFrom = xTo;
    }

    return result;
}