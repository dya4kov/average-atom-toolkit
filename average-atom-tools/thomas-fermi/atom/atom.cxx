#include <cmath>
#include <numeric>
#include <algorithm>

#include <numeric-tools/specfunc/fermi-dirac/complete.h>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/atom.h>
#include <average-atom-tools/thomas-fermi/atom/ODE/potential.h>

using namespace AATools::TF;
using numtools::ODE::Array;
using numtools::ODE::Solver;
using numtools::ODE::stepper::PD853;
using numtools::specfunc::FermiDirac;
using numtools::specfunc::FD::Half;
using AATools::TF::ODE::RHSPotential;

Atom::Atom() {
    V1 = 1.0; T1 = 1.0;
    VZ = 1.0; TZ = 1.0;
    tolerance = 1e-6;
    e .setTolerance(tolerance);
    RP.setTolerance(tolerance);
    xU.setTolerance(tolerance);
    muZ = xU.mu();
    mu1 = muZ;
}

void Atom::setTolerance(const double& t) {
    double Z = V1/VZ;
    tolerance = t;
    e .setTolerance(tolerance);
    RP.setTolerance(tolerance);
    xU.setTolerance(tolerance);
    muZ = xU.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

void Atom::setV(const double& V) {
    double Z = V1/VZ;
    VZ = V;
    V1 = V*Z;
    e .setV(V);
    RP.setV(V);
    xU.setV(V);
    muZ = xU.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

void Atom::setT(const double& T) {
    double Z = V1/VZ;
    TZ = T;
    T1 = T*std::pow(Z, -4.0/3.0);
    e .setT(T);
    RP.setT(T);
    xU.setT(T);
    muZ = xU.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

void Atom::setZ(const double& Z) {
    double Zold = V1/VZ;
    V1 = V1*Z/Zold;
    VZ = V1/Z;
    T1 = T1*std::pow(Z/Zold, -4.0/3.0);
    TZ = T1*std::pow(Z, 4.0/3.0);
    e .setZ(Z);
    RP.setZ(Z);
    xU.setZ(Z);
    muZ = xU.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
}

// double Atom::tauEval(const double& e, const double& l) {
// 
//     ::numtools::ODE::Array<RHSTau::dim> y; y.fill(0.0);
//     double rp2 = rPoint2(e, l);
//     double rp1 = rPoint1(e, l);
//     double result = 0.0;
//     double Z = V1/VZ;
// 
//     if (rp2 > rp1) {
//         y[0] = xU(rp2);
//         y[1] = xUDX(rp2);
// 
//         rhsTau.set_e(e);
//         rhsTau.set_l(l);
//         
//         tauSolver.setStep(1e-10);
//         tauSolver.integrate(rhsTau, y, std::sqrt(rp2), std::sqrt(rp1));
// 
//         result = 2.0*std::sqrt(2.0)*r0()*y[2]/Z;
//     }
// 
//     return result;
// }

// double Atom::tau(const int& n, const int& l) {
//     if (needAdjust) adjust();
// 
//     if (tLevel.size() < idxLevel(n + 1) + 1) {
//         tLevel.resize(idxLevel(n + 1) + 1);
//         tReady.resize(idxLevel(n + 1) + 1, false);
//     }
// 
//     double Z = V1/VZ;
// 
//     if (tReady[idxLevel(n) + l]) 
//         return tLevel[idxLevel(n) + l]*r0()*std::pow(Z, -1.0/3.0);
// 
//     double lambda = l + 0.5;
//     double lArg = 0.5*lambda*lambda / r0() / r0() * std::pow(Z, -2.0/3.0);
//     double eArg = e(n,l) * std::pow(Z, -4.0/3.0);
//     tLevel[idxLevel(n) + l] = tauEval(eArg, lArg);
// 
//     return tLevel[idxLevel(n) + l];//*r0()*std::pow(Z, -1.0/3.0);
// }

double Atom::eDens(const double& x) {
    
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

double* Atom::eDens(double* x, size_t n) {
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](size_t i1, size_t i2) {
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

std::vector<double>& Atom::eDens(const std::vector<double>& x) {
    size_t n = x.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double>* result = new std::vector<double>(n);

    std::sort(idx.begin(), idx.end(),
       [&x](size_t i1, size_t i2) {
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

double Atom::rpInner(const int& n, const int& l) {
    double Z = V1/VZ;
    double lambda = l + 0.5;
    double r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double eArg = std::pow(Z, -4.0/3.0)*e(n, l);
    return RP.inner(eArg, lArg);
}

double Atom::rpOuter(const int& n, const int& l) {
    double Z = V1/VZ;
    double lambda = l + 0.5;
    double r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double eArg = std::pow(Z, -4.0/3.0)*e(n, l);
    return RP.outer(eArg, lArg);
}