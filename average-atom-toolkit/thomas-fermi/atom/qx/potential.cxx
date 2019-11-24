#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/qx/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/qx/ODE/potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/qx/chemical-potential.h>

using namespace aatk::TF::qx;
using ::numtk::ODE::Array;
using ::numtk::ODE::Solver;
using ::numtk::ODE::stepper::PD853;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::MHalf;

using ::aatk::TF::qx::ODE::RHSPotential;

Potential::Potential() : 
    V(1.0), T(1.0), Z(1.0), 
    tolerance(1e-6)
{
    ::aatk::TF::ChemicalPotential M;
    ::aatk::TF::qx::ChemicalPotential dM;
    M.setZ(Z); dM.setZ(Z);
    M.setTolerance(tolerance);
    dM.setTolerance(tolerance);
    mu = M(V, T);
    dmuqx = dM(V, T);
}

Potential::Potential(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    mu = p.mu; dmuqx = p.dmuqx;
    tolerance  = p.tolerance;
}

Potential& Potential::operator=(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    mu = p.mu; dmuqx = p.dmuqx;
    tolerance  = p.tolerance;
    return *this;
}

void Potential::setV(const double& _V) { 
    if (std::abs(V - _V) > tolerance) {
        V = _V;
        ::aatk::TF::ChemicalPotential M;
        ::aatk::TF::qx::ChemicalPotential dM;
        M.setZ(Z); dM.setZ(Z);
        M.setTolerance(tolerance);
        dM.setTolerance(tolerance);
        mu = M(V, T);
        dmuqx = dM(V, T);
    }
}

void Potential::setT(const double& _T) {
    if (std::abs(T - _T) > tolerance) {
        T = _T;
        ::aatk::TF::ChemicalPotential M;
        ::aatk::TF::qx::ChemicalPotential dM;
        M.setZ(Z); dM.setZ(Z);
        M.setTolerance(tolerance);
        dM.setTolerance(tolerance);
        mu = M(V, T);
        dmuqx = dM(V, T);
    }
}

void Potential::setZ(const double& _Z) {
    if (std::abs(Z - _Z) > tolerance) {
        Z = _Z;
        ::aatk::TF::ChemicalPotential M;
        ::aatk::TF::qx::ChemicalPotential dM;
        M.setZ(Z); dM.setZ(Z);
        M.setTolerance(tolerance);
        dM.setTolerance(tolerance);
        mu = M(V, T);
        dmuqx = dM(V, T);
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
        ::aatk::TF::ChemicalPotential M;
        ::aatk::TF::qx::ChemicalPotential dM;
        M.setZ(Z); dM.setZ(Z);
        M.setTolerance(tolerance);
        dM.setTolerance(tolerance);
        mu = M(V, T);
        dmuqx = dM(V, T);
    }
}

void Potential::setTolerance(const double& t) {
    if (std::abs(tolerance - t) > 1e-15) {
        tolerance = std::max(t, tolerance);
        ::aatk::TF::ChemicalPotential M;
        ::aatk::TF::qx::ChemicalPotential dM;
        M.setZ(Z); dM.setZ(Z);
        M.setTolerance(tolerance);
        dM.setTolerance(tolerance);
        mu = M(V, T);
        dmuqx = dM(V, T);
    }
}

double Potential::operator()(const double& x) {
    RHSPotential rhs;
    FermiDirac<MHalf> FDmhalf;

    double    mu1 = mu*std::pow(Z, -4.0/3.0);
    double     V1 = V*Z;
    double     T1 = T*std::pow(Z, -4.0/3.0);
    double dmuqx1 = dmuqx*std::pow(Z, -2.0/3.0);
    double chi1;

    if (T1 <= 1e-10) chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - std::sqrt(mu1);
    else chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - 0.5*std::sqrt(T1)*FDmhalf(mu1/T1);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> chi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    chi.fill(0.0);
    chi[2] = chi[3] = chi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    if (xFrom > xTo) solver.integrate(rhs, chi, xFrom, xTo);

    double psi;
    if (T1 <= 1e-10) psi = std::sqrt(2.0)/(6.0*M_PI)*(chi[2] + std::sqrt(x*(chi[0] + x*mu1)));
    else psi = std::sqrt(2.0)/(6.0*M_PI)*(chi[2] + 0.5*x*std::sqrt(T1)*FDmhalf((chi[0]/(x*T1) + mu1/T1)));

    return (psi - x*dmuqx1)*std::pow(Z, 2.0/3.0);
}

double Potential::dx(const double& x) {
    RHSPotential rhs;
    FermiDirac<MHalf> FDmhalf;

    double    mu1 = mu*std::pow(Z, -4.0/3.0);
    double     V1 = V*Z;
    double     T1 = T*std::pow(Z, -4.0/3.0);
    double dmuqx1 = dmuqx*std::pow(Z, -2.0/3.0);
    double chi1;

    if (T1 <= 1e-10) chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - std::sqrt(mu1);
    else chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - 0.5*std::sqrt(T1)*FDmhalf(mu1/T1);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> chi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    chi.fill(0.0);
    chi[2] = chi[3] = chi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    if (xFrom > xTo) solver.integrate(rhs, chi, xFrom, xTo);
    return chi[3];
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
    FermiDirac<MHalf> FDmhalf;

    double    mu1 = mu*std::pow(Z, -4.0/3.0);
    double     V1 = V*Z;
    double     T1 = T*std::pow(Z, -4.0/3.0);
    double dmuqx1 = dmuqx*std::pow(Z, -2.0/3.0);
    double chi1;

    if (T1 <= 1e-10) chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - std::sqrt(mu1);
    else chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - 0.5*std::sqrt(T1)*FDmhalf(mu1/T1);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> chi;
    double xFrom = 1.0;
    chi.fill(0.0);
    chi[2] = chi[3] = chi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    solver.setStep(tolerance);

    for (auto i : idx) {
        // solver.integrate(rhs, chi, xFrom, xTo);
        // double xTo = std::sqrt(x[i]);
        // double psi;
        // if (T1 <= 1e-10) psi = std::sqrt(2.0)/(6.0*M_PI)*(chi[2] + std::sqrt(x[i]*(chi[0] + x[i]*mu1)));
        // else psi = std::sqrt(2.0)/(6.0*M_PI)*(chi[2] + 0.5*x[i]*std::sqrt(T1)*FDmhalf((chi[0]/(x[i]*T1) + mu1/T1)));
        // result[i] = (psi - x[i]*dmuqx1)*std::pow(Z, 2.0/3.0);
        result[i] = operator()(x[i]);
        // std::cout << "x = " << x[i] << ", phi = " << chi[0] << ", chi = " << chi[2] << ", psi = " << psi << ", xdU = " << result[i] << std::endl;
        // xFrom = xTo;
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
    FermiDirac<MHalf> FDmhalf;

    double    mu1 = mu*std::pow(Z, -4.0/3.0);
    double     V1 = V*Z;
    double     T1 = T*std::pow(Z, -4.0/3.0);
    double dmuqx1 = dmuqx1*std::pow(Z, -2.0/3.0);
    double chi1;

    if (T1 <= 1e-10) chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - std::sqrt(mu1);
    else chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - 0.5*std::sqrt(T1)*FDmhalf(mu1/T1);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> chi;
    double xFrom = 1.0;
    chi.fill(0.0);
    chi[2] = chi[3] = chi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, chi, xFrom, xTo);
        result[i] = chi[3];
        xFrom = xTo;
    }

    return result;
}

std::vector<double> Potential::operator()(const std::vector<double>& x) {

    std::size_t n = x.size();
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double> result(n);

    std::sort(idx.begin(), idx.end(),
       [&x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;
    FermiDirac<MHalf> FDmhalf;

    double    mu1 = mu*std::pow(Z, -4.0/3.0);
    double     V1 = V*Z;
    double     T1 = T*std::pow(Z, -4.0/3.0);
    double dmuqx1 = dmuqx1*std::pow(Z, -2.0/3.0);
    double chi1;

    if (T1 <= 1e-10) chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - std::sqrt(mu1);
    else chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - 0.5*std::sqrt(T1)*FDmhalf(mu1/T1);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> chi;
    double xFrom = 1.0;
    chi.fill(0.0);
    chi[2] = chi[3] = chi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, chi, xFrom, xTo);
        double psi;
        if (T1 <= 1e-10) psi = std::sqrt(2.0)/(6.0*M_PI)*(chi[2] + std::sqrt(x[i]*(chi[0] + x[i]*mu1)));
        else psi = std::sqrt(2.0)/(6.0*M_PI)*(chi[2] + 0.5*x[i]*std::sqrt(T1)*FDmhalf((chi[0] + x[i]*mu1)/(x[i]*T1)));
        result[i] = (psi/x[i] - dmuqx1)*std::pow(Z, 2.0/3.0);
        xFrom = xTo;
    }

    return result;
}

std::vector<double> Potential::dx(const std::vector<double>& x) {

    std::size_t n = x.size();
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double> result(n);

    std::sort(idx.begin(), idx.end(),
       [&x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;
    FermiDirac<MHalf> FDmhalf;

    double    mu1 = mu*std::pow(Z, -4.0/3.0);
    double     V1 = V*Z;
    double     T1 = T*std::pow(Z, -4.0/3.0);
    double dmuqx1 = dmuqx1*std::pow(Z, -2.0/3.0);
    double chi1;

    if (T1 <= 1e-10) chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - std::sqrt(mu1);
    else chi1 = (6.0 * M_PI)/std::sqrt(2.0)*dmuqx1 - 0.5*std::sqrt(T1)*FDmhalf(mu1/T1);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSPotential::dim> chi;
    double xFrom = 1.0;
    chi.fill(0.0);
    chi[2] = chi[3] = chi1;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, chi, xFrom, xTo);
        result[i] = chi[3];
        xFrom = xTo;
    }

    return result;
}