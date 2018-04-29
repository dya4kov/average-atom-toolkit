#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/qx/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/qx/ODE/potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/qx/chemical-potential.h>

using namespace aatk::TF::qx;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using aatk::TF::qx::ODE::RHSPotential;

Potential::Potential() : 
    V(1.0), T(1.0), Z(1.0),
    psi1(1.0),
    tolerance(1e-6),
    bestTolerance(1e-12),
    needUpdate(true) 
{
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
}

Potential::Potential(const Potential& p) : bestTolerance(1e-12) {
    V = p.V; T = p.T; Z = p.Z;
    psi1 = p.psi1;
    tolerance  = p.tolerance;
    needUpdate = p.needUpdate;
    mu = p.mu;
}

Potential& Potential::operator=(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    psi1 = p.psi1;
    tolerance  = p.tolerance;
    needUpdate = p.needUpdate;
    mu = p.mu;
    return *this;
}

void Potential::setV(const double& _V) { 
    if (std::abs(V - _V) > bestTolerance) {
        V = _V;
        ChemicalPotential M;
        M.setZ(Z);
        M.setTolerance(tolerance);
        mu = M(V, T);
        needUpdate = true;
    }
}

void Potential::setT(const double& _T) {
    if (std::abs(T - _T) > bestTolerance) {
        T = _T;
        ChemicalPotential M;
        M.setZ(Z);
        M.setTolerance(tolerance);
        mu = M(V, T);
        needUpdate = true;
    }
}

void Potential::setZ(const double& _Z) {
    if (std::abs(Z - _Z) > bestTolerance) {
        Z = _Z;
        ChemicalPotential M;
        M.setZ(Z);
        M.setTolerance(tolerance);
        mu = M(V, T);
        needUpdate = true;
    }
}

void Potential::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    if (std::abs(V - _V) > bestTolerance ||
        std::abs(T - _T) > bestTolerance ||
        std::abs(Z - _Z) > bestTolerance) {
        V = _V; T = _T; Z = _Z;
        ChemicalPotential M;
        M.setZ(Z);
        M.setTolerance(tolerance);
        mu = M(V, T);
        needUpdate = true;
    }
}

void Potential::setTolerance(const double& t) {
    if (std::abs(tolerance - t) > bestTolerance) {
        tolerance = std::max(t, bestTolerance);
        ChemicalPotential M;
        M.setZ(Z);
        M.setTolerance(tolerance);
        mu = M(V, T);
        needUpdate = true;
    }
}

double Potential::operator()(const double& x) {
    RHSPotential rhs;

    if (needUpdate) update();

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    psi.fill(0.0);
    psi[2] = psi[3] = psi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    if (xFrom > xTo) solver.integrate(rhs, psi, xFrom, xTo);
    return psi[2];
}

double Potential::dx(const double& x) {
    RHSPotential rhs;

    if (needUpdate) update();

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    psi.fill(0.0);
    psi[2] = psi[3] = psi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    if (xFrom > xTo) solver.integrate(rhs, psi, xFrom, xTo);
    return psi[3];
}

double* Potential::operator()(const double* x, const std::size_t& n) {
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;

    if (needUpdate) update();

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> psi; 
    double xFrom = 1.0;
    psi.fill(0.0);
    psi[2] = psi[3] = psi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        result[i] = psi[2];
        xFrom = xTo;
    }

    return result;
}

double* Potential::dx(const double* x, const std::size_t& n) {
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;

    if (needUpdate) update();

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> psi; 
    double xFrom = 1.0;
    psi.fill(0.0);
    psi[2] = psi[3] = psi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        result[i] = psi[3];
        xFrom = xTo;
    }

    return result;
}

std::vector<double>& Potential::operator()(const std::vector<double>& x) {

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

    if (needUpdate) update();

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;
    psi.fill(0.0);
    psi[2] = psi[3] = psi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        (*result)[i] = psi[2];
        xFrom = xTo;
    }

    return *result;
}

std::vector<double>& Potential::dx(const std::vector<double>& x) {

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

    if (needUpdate) update();

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> psi; 
    double xFrom = 1.0;
    psi.fill(0.0);
    psi[2] = psi[3] = psi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, psi, xFrom, xTo);
        (*result)[i] = psi[3];
        xFrom = xTo;
    }

    return *result;
}