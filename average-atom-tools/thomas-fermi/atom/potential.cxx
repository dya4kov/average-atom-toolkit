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
    V1(1.0), T1(1.0), mu1(4.10057773),
    VZ(1.0), TZ(1.0), muZ(4.10057773), 
    tolerance(1e-6), 
    bestTolerance(1e-12),
    needUpdate(false) {}

Potential::Potential(const Potential& p) : bestTolerance(1e-12) {
    V1 = p.V1; T1 = p.T1; mu1 = p.mu1;
    VZ = p.VZ; TZ = p.TZ; muZ = p.muZ;
    tolerance  = p.tolerance;
    needUpdate = p.needUpdate;
}

Potential& Potential::operator=(const Potential& p) {
    V1 = p.V1; T1 = p.T1; mu1 = p.mu1;
    VZ = p.VZ; TZ = p.TZ; muZ = p.muZ;
    tolerance  = p.tolerance;
    needUpdate = p.needUpdate;
    return *this;
}

void Potential::setV(const double& V) { 
    if (std::abs(VZ - V) > bestTolerance) {
        double Z = V1/VZ;
        VZ = V; V1 = V*Z; 
        needUpdate = true;
    }
}

void Potential::setT(const double& T) {
    if (std::abs(TZ - T) > bestTolerance) {
        double Z = V1/VZ;
        TZ = T; T1 = T*std::pow(Z, -4.0/3.0); 
        needUpdate = true;
    }
}

void Potential::setZ(const double& Z) {
    double Zold = V1/VZ;
    if (std::abs(Z - Zold) > bestTolerance) {
        T1  = TZ*std::pow(Z, -4.0/3.0);
        V1  = VZ*Z;
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

    if (needUpdate) mu();

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
    RHSPotential rhs;

    if (needUpdate) mu();

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

    RHSPotential rhs;

    if (needUpdate) mu();

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

    RHSPotential rhs;

    if (needUpdate) mu();

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

    RHSPotential rhs;

    if (needUpdate) mu();

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

    RHSPotential rhs;

    if (needUpdate) mu();

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

double Potential::mu() {
	if (!needUpdate) return muZ;

    Array<RHSPotential::dim> phi; 

    RHSPotential rhs;
	rhs.set_V(V1);
    rhs.set_T(T1);
    
    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    double phi_0 = std::pow(4.0*M_PI/3.0/V1, 1.0/3.0);

    double xFrom = 1.0;
    double xTo   = 0.0;
    double delta = 0.1;

    if (T1 <= 1e-10) T1 = 1e-11;
    double muStart = table(std::log10(V1), std::log10(T1));

    bool convergenceSuccess = false;
    while (!convergenceSuccess) {
        double muLeftStart  = muStart - delta*std::max(10.0, std::abs(muStart));
        double muRightStart = muStart + delta*std::max(10.0, std::abs(muStart));
        double muLeft  = muLeftStart;
        double muRight = muRightStart;
        double error   = std::abs(muRight - muLeft)/std::abs(muLeft + muRight + tolerance);
        while (error > tolerance) {
            phi.fill(0.0);
            rhs.set_mu(0.5*(muLeft + muRight));

            solver.integrate(rhs, phi, xFrom, xTo);

            if (std::isfinite(phi[0])) {
               if (phi[0] - phi_0 > 0)
                   muRight -= 0.5*(muRight - muLeft);
               else muLeft += 0.5*(muRight - muLeft);
            }
            else {
                muRight -= 0.5*(muRight - muLeft);
            }

            error = std::abs(muRight - muLeft)/std::abs(muRight + muLeft);
        }
        mu1 = 0.5*(muLeft + muRight);
        convergenceSuccess = true;
        if (std::abs( muLeftStart - muLeft )/std::abs( muLeftStart + muLeft ) < tolerance*100.0 || 
            std::abs(muRightStart - muRight)/std::abs(muRightStart + muRight) < tolerance*100.0  ) {
            convergenceSuccess = false;
            delta = delta + 0.1;
            if (delta > 0.45) {
                std::cout << "Potential->mu() error: too big delta" << std::endl;
                exit(0);
            }
        }
    }
    double Z = V1/VZ;
    muZ = mu1*std::pow(Z, 4.0/3.0);
    needUpdate = false;
    return muZ;
}