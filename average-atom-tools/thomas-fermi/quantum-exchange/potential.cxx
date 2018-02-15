#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/quantum-exchange/potential.h>
#include <average-atom-tools/thomas-fermi/quantum-exchange/ODE/potential.h>

using namespace AATools::TF::QE;
using numtools::ODE::Array;
using numtools::ODE::Solver;
using numtools::ODE::stepper::PD853;

using AATools::TF::QE::ODE::RHSPotential;

Potential::Potential() : 
    V1(1.0), T1(1.0), mu1(4.10057773),
    VZ(1.0), TZ(1.0), muZ(4.10057773),
    psi1(1.0),
    tolerance(1e-6),
    bestTolerance(1e-12),
    needUpdate(true) {}

Potential::Potential(const Potential& p) : bestTolerance(1e-12) {
    V1 = p.V1; T1 = p.T1; mu1 = p.mu1;
    VZ = p.VZ; TZ = p.TZ; muZ = p.muZ;
    psi1 = p.psi1;
    tolerance  = p.tolerance;
    needUpdate = p.needUpdate;
    phi = p.phi;
}

Potential& Potential::operator=(const Potential& p) {
    V1 = p.V1; T1 = p.T1; mu1 = p.mu1;
    VZ = p.VZ; TZ = p.TZ; muZ = p.muZ;
    psi1 = p.psi1;
    tolerance  = p.tolerance;
    needUpdate = p.needUpdate;
    phi = p.phi;
    return *this;
}

void Potential::setV(const double& V) { 
    if (std::abs(VZ - V) > bestTolerance) {
        double Z = V1/VZ;
        VZ = V; V1 = V*Z;
        phi.setV(V); 
        muZ = phi.mu();
        mu1 = muZ*std::pow(Z, -4.0/3.0);
        needUpdate = true;
    }
}

void Potential::setT(const double& T) {
    if (std::abs(TZ - T) > bestTolerance) {
        double Z = V1/VZ;
        TZ = T; T1 = T*std::pow(Z, -4.0/3.0);
        phi.setT(T);
        muZ = phi.mu();
        mu1 = muZ*std::pow(Z, -4.0/3.0);
        needUpdate = true;
    }
}

void Potential::setZ(const double& Z) {
    double Zold = V1/VZ;
    if (std::abs(Z - Zold) > bestTolerance) {
        T1  = TZ*std::pow(Z, -4.0/3.0);
        V1  = VZ*Z;
        phi.setZ(Z);
        muZ = phi.mu();
        mu1 = muZ*std::pow(Z, -4.0/3.0);
        needUpdate = true;
    }
}

void Potential::setTolerance(const double& t) {
    if (std::abs(tolerance - t) > bestTolerance) {
        tolerance = std::max(t, bestTolerance);
        needUpdate = true;
    }
}

double Potential::operator()(const double& x) {
    RHSPotential rhs;

    if (needUpdate) update();

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

double* Potential::operator()(const double* x, const size_t& n) {
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](size_t i1, size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;

    if (needUpdate) update();

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

double* Potential::dx(const double* x, const size_t& n) {
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](size_t i1, size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;

    if (needUpdate) update();

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

    if (needUpdate) update();

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

    if (needUpdate) update();

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

void Potential::update() {

    if (!needUpdate) return;
    
    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    Array<RHSPotential::dim> psi;
    psi.fill(0.0);

    RHSPotential rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    double xFrom = 1.0;
    double xTo   = 0.0;
    double delta = 0.1;
    bool convergenceSuccess = false;

    if (T1 <= 1e-10) T1 = 1e-11;
    double psiStart = table(std::log10(V1), std::log10(T1));
    while (!convergenceSuccess) {
        double psiLeftStart  = psiStart - delta*std::abs(psiStart);
        double psiRightStart = psiStart + delta*std::abs(psiStart);
        double psiLeft  = psiLeftStart;
        double psiRight = psiRightStart;
        double error   = std::abs(psiRight - psiLeft)/std::abs(psiLeft + psiRight + tolerance);
        while (error > tolerance) {
            psi[0] = 0.0;
            psi[1] = 0.0;
            psi[2] = 0.5*(psiLeft + psiRight);
            psi[3] = psi[2];

            solver.integrate(rhs, psi, xFrom, xTo);

            if (psi[2] > 0)
                psiRight -= 0.5*(psiRight - psiLeft);
            else psiLeft += 0.5*(psiRight - psiLeft);

            error = std::abs(psiRight - psiLeft)/std::abs(psiLeft + psiRight);
        }
        psi1 = 0.5*(psiLeft + psiRight);
        convergenceSuccess = true;
        if (std::abs( psiLeftStart - psiLeft )/std::abs( psiLeftStart + psiLeft ) < tolerance*100.0 || 
            std::abs(psiRightStart - psiRight)/std::abs(psiRightStart + psiRight) < tolerance*100.0  ) {
            convergenceSuccess = false;
            delta = delta + 0.1;
            if (delta > 0.45) {
                std::cout << "potential error: too big delta" << std::endl;
                exit(0);
            }
        }
    }
    needUpdate = false;
}