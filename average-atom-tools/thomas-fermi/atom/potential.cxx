#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/atom/potential.h>
#include <average-atom-tools/thomas-fermi/atom/ODE/potential.h>

using namespace AATools::TF;
using numtools::ODE::Array;
using numtools::ODE::Solver;
using numtools::ODE::stepper::PD853;
using AATools::TF::ODE::RHSPotential;

Potential::Potential() : 
    V(1.0), T(1.0), Z(1.0),
    tolerance(1e-6) {}

Potential::Potential(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    tolerance  = p.tolerance;
}

Potential& Potential::operator=(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    tolerance  = p.tolerance;
    return *this;
}

void Potential::setV(const double& _V) { 
    V = _V;
}

void Potential::setT(const double& _T) {
    T = _T;
}

void Potential::setZ(const double& _Z) {
    Z = _Z; mu.setZ(Z);
}

void Potential::setTolerance(const double& t) {
    tolerance = t; mu.setTolerance(t);
}

double Potential::operator()(const double& x) {

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

    solver.integrate(rhs, phi, xFrom, xTo);
    return phi[0];
}

double Potential::dx(const double& x) {

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
    solver.integrate(rhs, phi, xFrom, xTo);
    return phi[1];
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

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, phi, xFrom, xTo);
        result[i] = phi[0];
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

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, phi, xFrom, xTo);
        result[i] = phi[1];
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

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, phi, xFrom, xTo);
        (*result)[i] = phi[0];
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

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, phi, xFrom, xTo);
        (*result)[i] = phi[1];
        xFrom = xTo;
    }

    return *result;
}
