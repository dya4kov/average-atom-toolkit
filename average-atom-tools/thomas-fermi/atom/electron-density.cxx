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

ElectronDensity::ElectronDensity() :
    V(1.0), T(1.0), Z(1.0),
    tolerance(1e-6)
{
    mu.setTolerance(tolerance);
}

void ElectronDensity::setTolerance(const double& t) {
    mu.setTolerance(tolerance);
}

void ElectronDensity::setV(const double& _V) {
    V = _V;
}

void ElectronDensity::setT(const double& _T) {
    T = _T;
}

void ElectronDensity::setZ(const double& _Z) {
    Z = _Z; mu.setZ(Z);
}

double ElectronDensity::operator()(const double& x) {

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);
    
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
    return T*std::sqrt(2.0*T)*FDhalf(phi[0]/(x*T1) + mu1/T1)/(M_PI*M_PI);
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

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

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
        result[i] = T*std::sqrt(2.0*T)*FDhalf(phi[0]/(x[i]*T1) + mu1/T1)/(M_PI*M_PI);
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

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

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
        (*result)[i] = T*std::sqrt(2.0*T)*FDhalf(phi[0]/(x[i]*T1) + mu1/T1)/(M_PI*M_PI);
        xFrom = xTo;
    }

    return *result;
}