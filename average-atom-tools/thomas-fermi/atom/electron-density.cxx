#include <cmath>
#include <numeric>
#include <algorithm>

#include <numeric-tools/specfunc/fermi-dirac/complete.h>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/atom/electron-density.h>
#include <average-atom-tools/thomas-fermi/atom/ODE/potential.h>

using namespace AATools::TF;
using numtools::ODE::Array;
using numtools::ODE::Solver;
using numtools::ODE::stepper::PD853;
using numtools::specfunc::FermiDirac;
using numtools::specfunc::FD::Half;
using AATools::TF::ODE::RHSPotential;

ElectronDensity::ElectronDensity() {
    V1 = 1.0; T1 = 1.0;
    VZ = 1.0; TZ = 1.0;
    tolerance = 1e-6;
    potential.setTolerance(tolerance);
    muZ = potential.mu();
    mu1 = muZ;
}

void ElectronDensity::setTolerance(const double& t) {
    double Z = V1/VZ;
    tolerance = t;
    potential.setTolerance(tolerance);
    muZ = potential.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

void ElectronDensity::setV(const double& V) {
    double Z = V1/VZ;
    VZ = V;
    V1 = V*Z;
    potential.setV(V);
    muZ = potential.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

void ElectronDensity::setT(const double& T) {
    double Z = V1/VZ;
    TZ = T;
    T1 = T*std::pow(Z, -4.0/3.0);
    potential.setT(T);
    muZ = potential.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

void ElectronDensity::setZ(const double& Z) {
    double Zold = V1/VZ;
    V1 = V1*Z/Zold;
    VZ = V1/Z;
    T1 = T1*std::pow(Z/Zold, -4.0/3.0);
    TZ = T1*std::pow(Z, 4.0/3.0);
    potential.setZ(Z);
    muZ = potential.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

double ElectronDensity::operator()(const double& x) {
    
    RHSPotential rhs;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> phi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    phi.fill(0.0);

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    FermiDirac<Half> FDhalf;

    if (xTo < xFrom) solver.integrate(rhs, phi, xFrom, xTo);
    return TZ*std::sqrt(2.0*TZ)*FDhalf(phi[0]/(x*T1) + mu1/T1)/(M_PI*M_PI);
}

double* ElectronDensity::operator()(const double* x, const std::size_t& n) {
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> phi; 
    double xFrom = 1.0;
    phi.fill(0.0);

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    FermiDirac<Half> FDhalf;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, phi, xFrom, xTo);
        result[i] = TZ*std::sqrt(2.0*TZ)*FDhalf(phi[0]/(x[i]*T1) + mu1/T1)/(M_PI*M_PI);
        xFrom = xTo;
    }

    return result;
}

std::vector<double>& ElectronDensity::operator()(const std::vector<double>& x) {
    std::size_t n = x.size();
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double>* result = new std::vector<double>(n);

    std::sort(idx.begin(), idx.end(),
       [&x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> phi; 
    double xFrom = 1.0;
    phi.fill(0.0);

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    FermiDirac<Half> FDhalf;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, phi, xFrom, xTo);
        (*result)[i] = TZ*std::sqrt(2.0*TZ)*FDhalf(phi[0]/(x[i]*T1) + mu1/T1)/(M_PI*M_PI);
        xFrom = xTo;
    }

    return *result;
}