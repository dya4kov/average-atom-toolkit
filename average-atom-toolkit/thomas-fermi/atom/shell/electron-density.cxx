#include <cmath>
#include <numeric>
#include <algorithm>

//#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/wave-function.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/electron-density.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

using namespace aatk::TF::shell;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

ElectronDensity::ElectronDensity() :
    V(1.0), T(1.0), Z(1.0), mu(4.100577730112),
    tolerance(1e-6), eBoundary(1e+15), 
    nmax(15), ready(false) {}

void ElectronDensity::setTolerance(const double& t) {
    tolerance = t;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setTolerance(t);
    ready = false;
}

void ElectronDensity::setThreadsLimit(const int& Nthreads) {
    e.setThreadsLimit(Nthreads);
}
void ElectronDensity::setNmax(const int& Nmax) { nmax = Nmax; }

void ElectronDensity::setBoundary(const double& eb) { eBoundary = eb; }

void ElectronDensity::setV(const double& _V) {
    V = _V;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setV(V);
    ready = false;
}

void ElectronDensity::setT(const double& _T) {
    T = _T;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setT(T);
    ready = false;
}

void ElectronDensity::setZ(const double& _Z) {
    Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setZ(Z);
    ready = false;
}

void ElectronDensity::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    V = _V; T = _T; Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setVTZ(V, T, Z);
    ready = false;
}

double ElectronDensity::operator()(const double& x) {

    if (!ready) e.prepareLevelsBelow(nmax);
    ready = true;

    double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);

    ::aatk::TF::shell::WaveFunction WF;
    WF.setVTZ(V, T, Z);
    WF.setTolerance(tolerance);

    double rho = 0.0;

    for (int n = 1; n < nmax; ++n) {
        for (int l = 0; l < n; ++l) {
            double enl    = e(n,l);
            if (enl > eBoundary) continue;
            double Nnl    = (2.0*l + 1)/(1.0 + std::exp((enl - mu)/T));
            double lambda = l + 0.5;
            double Rnl    = WF(enl, lambda, x);
            rho          += Nnl*Rnl*Rnl;
        }
    }

    r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
    rho *= 1.0/(2.0*M_PI*r0*r0*x*x);

    return rho;
}

double* ElectronDensity::operator()(const double* x, const std::size_t& size) {

    if (!ready) e.prepareLevelsBelow(nmax);
    ready = true;

    double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);

    ::aatk::TF::shell::WaveFunction Rnl;
    Rnl.setVTZ(V, T, Z);
    Rnl.setTolerance(tolerance);

    double* rho = new double[size];
    for (std::size_t i = 0; i < size; ++i) rho[i] = 0.0;

    for (int n = 1; n < nmax; ++n) {
        for (int l = 0; l < n; ++l) {
            double  enl    = e(n,l);
            if (enl > eBoundary) continue;
            double  Nnl    = (2.0*l + 1.0)/(1.0 + std::exp((enl - mu)/T));
            double  lambda = l + 0.5;
            double* wf     = Rnl(enl, lambda, x, size);
            for (std::size_t i = 0; i < size; ++i) {
                rho[i] += 1.0/(x[i]*x[i])*Nnl*wf[i]*wf[i];
            }
            delete[] wf;
        }
    }

    r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
    for (std::size_t i = 0; i < size; ++i) {
        rho[i] *= 1.0/(2.0*M_PI*r0*r0);
    }

    return rho;
}

std::vector<double> ElectronDensity::operator()(const std::vector<double>& x) {
    double*  data = operator()(x.data(), x.size());
    std::vector<double> result(data, data + x.size());
    delete[] data;
    return result;
}