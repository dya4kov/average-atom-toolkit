#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/potential-old.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/ODE/ksi.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/ODE/potential-old.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/shell/ODE/dM.h>

using namespace aatk::TF::shell;
using ::numtk::ODE::Array;
using ::numtk::ODE::Solver;
using ::numtk::ODE::stepper::PD853;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::MHalf;

using ::aatk::TF::shell::ODE::RHSksi;
using ::aatk::TF::shell::ODE::RHSPotential;
using ::aatk::TF::shell::ODE::RHSdM;

Potential::Potential() : 
    V(1.0), T(1.0), Z(1.0), nmax(25), 
    tolerance(1e-6)
{
    N.setNmax(nmax);
    ready = false;
}

Potential::Potential(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    mu = p.mu; dmush = p.dmush;
    N = p.N; nmax = p.nmax;
    eBoundary = p.eBoundary;
    tolerance = p.tolerance;
    ready = p.ready;
}

Potential& Potential::operator=(const Potential& p) {
    V = p.V; T = p.T; Z = p.Z;
    mu = p.mu; dmush = p.dmush;
    N = p.N; nmax = p.nmax;
    eBoundary = p.eBoundary;
    tolerance = p.tolerance;
    ready = p.ready;
    return *this;
}

void Potential::setV(const double& _V) { 
    if (std::abs(V - _V) > tolerance) {
        V = _V;
        N.setV(V); 
        ready = false;
    }
}

void Potential::setT(const double& _T) {
    if (std::abs(T - _T) > tolerance) {
        T = _T; 
        N.setT(T);
        ready = false;
    }
}

void Potential::setZ(const double& _Z) {
    if (std::abs(Z - _Z) > tolerance) {
        Z = _Z; 
        N.setZ(Z);
        ready = false;
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
        N.setVTZ(V, T, Z);
        ready = false;
    }
}

void Potential::setTolerance(const double& t) {
    if (std::abs(tolerance - t) > 1e-15) {
        tolerance = std::max(t, tolerance);
        N.setTolerance(t);
        ready = false;
    }
}

void Potential::setNmax(const int& _nmax) {
    if (_nmax != nmax) {
        nmax = _nmax; 
        N.setNmax(nmax);
        ready = false;
    }
}

void Potential::prepare() {

    // resize arrays and fill with 0
    e    .resize(RHSPotential::dim, 0.0);
    rpi  .resize(RHSPotential::dim, 0.0);
    rpo  .resize(RHSPotential::dim, 0.0);
    ksi0 .resize(RHSPotential::dim, 0.0);
    ksi21.resize(RHSPotential::dim, 0.0);
    sign .resize(RHSPotential::dim, 0.0);

    // prepare energy levels
    eBoundary  = N.eBoundary();

    double continuous = N.continuous(eBoundary);
    double discrete   = N.discrete(eBoundary);
    double dN = discrete - continuous;

    ::aatk::TF::ChemicalPotential M;
    M.setTolerance(tolerance);
    M.setZ(Z);

    auto eLevel = N.eLevel();

    mu = M(V, T);

    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);
    double mu1 = mu*std::pow(Z, -4.0/3.0);

    double  r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);

    // chemical potential
    
    RHSdM rhsdM;

    rhsdM.set_V(V1);
    rhsdM.set_T(T1);
    rhsdM.set_mu(mu1);

    Array<RHSdM::dim> ydM; ydM.fill(0.0);

    Solver<PD853<RHSdM>> solverdM;
    solverdM.setTolerance(0.0, 0.1*tolerance);

    solverdM.integrate(rhsdM, ydM, 1.0, 0.0);
    double dM = ydM[2]*6.0*V1*std::sqrt(2.0)/(M_PI*M_PI)*std::pow(Z, -1.0/3.0);
    
    dmush = -dN/dM;

    // rotate points, ksi0, ksi21, sign

    ::aatk::TF::RotatePoints RP;
    RP.setVTZ(V, T, Z);
    RP.setTolerance(tolerance);

    RHSksi rhsKsi;
    rhsKsi.set_V(V1);
    rhsKsi.set_T(T1);
    rhsKsi.set_mu(mu1);

    Solver<PD853<RHSksi>> solverKsi;
    solverKsi.setTolerance(0.1*tolerance, 0.0);

    ::numtk::specfunc::bessel::Jnu Jnu;
    ::numtk::specfunc::bessel::Ynu Ynu;
    ::numtk::specfunc::bessel::Knu Knu;

    for (int n = 1; n < nmax; ++n) {
        for (int l = 0; l < n; ++l) {

            int i = 2 + n*(n - 1)/2 + l;
            e[i] = eLevel(n, l);

            rpo[i] = std::sqrt(RP.outer(e[i], l + 0.5));
            rpi[i] = std::sqrt(RP.inner(e[i], l + 0.5));

            if (rpo[i] < tolerance) e[i] = 1e+15;

            if (e[i] > eBoundary) continue;

            double lambda = 0.5*(l + 0.5)*(l + 0.5) / r0 / r0 * std::pow(Z, -2.0/3.0);

            rhsKsi.set_e(e[i]*std::pow(Z, -4.0/3.0));
            rhsKsi.set_l(lambda);

            ksi0[i] = 0.0;

            Array<RHSksi::dim> yksi;

            if (rpo[i] < 1.0) {
                yksi[0] = 0.0; 
                yksi[1] = 0.0;
                yksi[2] = 0.0;
                
                solverKsi.setStep(tolerance);
                solverKsi.integrate(rhsKsi, yksi, 1.0, rpo[i]);

                ksi0[i] = 2.0*r0*std::sqrt(2.0)*yksi[2]*std::pow(Z, 1.0/3.0);

            }

            yksi[0] = RP.outerY(e[i], l + 0.5)[0];
            yksi[1] = RP.outerY(e[i], l + 0.5)[1];
            yksi[2] = 0.0;

            solverKsi.setStep(tolerance);
            solverKsi.integrate(rhsKsi, yksi, rpo[i], rpi[i]);

            ksi21[i]    = 1e-8 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2])*std::pow(Z, 1.0/3.0);
            double Jp13 = Jnu(1.0/3.0, ksi21[i]);
            double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21[i]));
            sign[i]     = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;

            e[i] *= std::pow(Z, -4.0/3.0);

            // std::cout << i << "e = " << e[i] << ", rpi = " << rpi[i] << ", rpo = " << rpo[i] << ", ksi0 = " << ksi0[i] << ", ksi21 = " << ksi21[i] << std::endl;
        }
    }
    eBoundary *= std::pow(Z, -4.0/3.0);
    ready = true;
}

double Potential::operator()(const double& x) {
    
    if (!ready) prepare();

    RHSPotential rhs;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V     (V1);
    rhs.set_T     (T1);
    rhs.set_Z     (Z);
    rhs.set_mu    (mu1);

    rhs.set_eb    (eBoundary);
    rhs.set_nmax  (nmax);
    rhs.set_e     (e);
    rhs.set_ksi0  (ksi0);
    rhs.set_ksi21 (ksi21);
    rhs.set_RP    (rpi, rpo);
    rhs.set_sign  (sign);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);

    psi.fill(0.0);
    psi[RHSPotential::dim - 1] = psi[RHSPotential::dim - 2] = dmush;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(tolerance, 0.0);

    if (xFrom > xTo) solver.integrate(rhs, psi, xFrom, xTo);

    return psi[RHSPotential::dim - 2]/x - dmush;
}

double Potential::dx(const double& x) {
    if (!ready) prepare();

    return 0.0;
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

    for (auto i : idx) { result[i] = 0.0; }

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

    for (auto i : idx) { result[i] = 0.0; } 

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

    for (auto i : idx) {
        result[i] = 0.0;
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

    for (auto i : idx) {
        result[i] = 0.0;
    }

    return result;
}