#include <cmath>
#include <algorithm>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/atom/electron-states.h>
#include <average-atom-tools/thomas-fermi/atom/ODE/boundary-energy.h>
#include <average-atom-tools/thomas-fermi/atom/ODE/continuous-states.h>

using namespace AATools::TF;

using numtools::ODE::Array;
using numtools::ODE::Solver;
using numtools::ODE::stepper::PD853;

using AATools::TF::ODE::RHSBE;
using AATools::TF::ODE::RHSCS;
using AATools::TF::ODE::RHSCSF;

ElectronStates::ElectronStates() : 
    V1(1.0), T1(1.0), mu1(4.10057773),
    VZ(1.0), TZ(1.0), muZ(4.10057773),
    muShift(0.0), nMax(15),
    tolerance(1e-6)
{
    phi.setTolerance(tolerance);
    e.setTolerance(tolerance);
}

void ElectronStates::setTolerance(const double& eps) {
    tolerance = eps;
    phi.setTolerance(eps);
    e.setTolerance(eps);
}

void ElectronStates::setV(const double& V) {
    double Z = V1/VZ;
    VZ = V;
    V1 = V*Z;
    phi.setV(V);
    muZ = phi.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
    e.setV(V); 
}

void ElectronStates::setT(const double& T) {
    double Z = V1/VZ;
    TZ = T;
    T1 = T*std::pow(Z, -4.0/3.0);
    phi.setT(T);
    muZ = phi.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
    e.setT(T);
}

void ElectronStates::setZ(const double& Z) {
    double Zold = V1/VZ;
    V1 = V1*Z/Zold;
    VZ = V1/Z;
    T1 = T1*std::pow(Z/Zold, -4.0/3.0);
    TZ = T1*std::pow(Z, 4.0/3.0);
    phi.setZ(Z);
    muZ = phi.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);
    e.setZ(Z);
}

void ElectronStates::setNmax(const int& N) {
    nMax = N;
}

void ElectronStates::setMuShift(const double& dmu) {
    double Z = V1/VZ;
    muShift = dmu*std::pow(Z, -4.0/3.0);
}

double ElectronStates::operator()(const int& n, const int& l) {
    double Z = V1/VZ;
    double Nnl = 0.0;
    double enl = e(n, l)*std::pow(Z, -4.0/3.0);
    if (T1 > 1e-10) 
         Nnl = (2.0 * l + 1.0) / (1.0 + std::exp((enl - mu1 - muShift)/T1));
    else Nnl = enl < mu1 + muShift ? 2.0*l + 1.0 : 0;
    return 2.0*Nnl;
}

double ElectronStates::operator()(const int& n) {
    double Z = V1/VZ;
    double Nn = 0.0;
    auto en = e[n];
    for (int l = 0; l < n; ++l) {
        double Nnl;
        double enl = en[l]*std::pow(Z, -4.0/3.0);
        if (T1 > 1e-10) 
             Nnl = (2.0 * l + 1.0) / (1.0 + std::exp((enl - mu1 - muShift)/T1));
        else Nnl = enl < mu1 + muShift ? 2.0*l + 1.0 : 0;
        Nn += Nnl;
    }
    return 2.0*Nn;
}

double ElectronStates::discrete() {
    double N = 0.0;
    for (int n = 1; n <= nMax; ++n)
        N += operator()(n);
    return N;
}

double ElectronStates::discrete(const double& energy) {
    double N = 0.0;
    double Z = V1/VZ;
    for (int n = 1; n <= nMax; ++n) {
        auto en = e[n];
        double Nn = 0.0;
        for (int l = 0; l < n; ++l) {
            double enl = en[l];
            if (enl < energy) {
                double Nnl;
                enl *= std::pow(Z, -4.0/3.0);
                if (T1 > 1e-10) 
                     Nnl = (2.0 * l + 1.0) / (1.0 + std::exp((enl - mu1 - muShift)/T1));
                else Nnl = en[l] < mu1 + muShift ? 2.0*l + 1.0 : 0;
                Nn += Nnl;
            }
        }
        N += Nn;
    }
    return 2.0*N;
}

std::vector<double>& ElectronStates::discrete(const std::vector<double>& energy) {

    auto N = new std::vector<double>(energy.size());

    for (std::size_t i = 0; i < energy.size(); ++i) {
        (*N)[i] = discrete(energy[i]);
    }

    return *N;
}

double* ElectronStates::discrete(const double* energy, const std::size_t& n) {

    double* N = new double[n];

    for (std::size_t i = 0; i < n; ++i) {
        N[i] = discrete(energy[i]);
    }

    return N;
}

double ElectronStates::continuous() {
    Array<RHSCSF::dim> y;
    y.fill(0.0);
    double Z = V1/VZ;

    RHSCSF rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    rhs.set_dmu(muShift);

    Solver<PD853<RHSCSF>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    solver.integrate(rhs, y, 1.0, 0.0);

    double states = 4.0 * y[RHSCSF::result] * std::sqrt(2.0) * V1 / (M_PI * M_PI ) * Z;

    if (T1 > 1e-10) states *= 1.5*T1*std::sqrt(T1);
    return states;
}

double ElectronStates::continuous(const double& energy) {
    Array<RHSCS::dim> y;
    y.fill(0.0);
    double Z = V1/VZ;

    RHSCS rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    rhs.set_dmu(muShift);
    rhs.set_e(energy*std::pow(Z, -4.0/3.0));

    Solver<PD853<RHSCS>> solver;
    solver.setTolerance(0.0, 10.0*tolerance);
    solver.integrate(rhs, y, 1.0, 0.0);

    double CS = 4.0 * y[RHSCS::result] * std::sqrt(2.0) * V1 / (M_PI * M_PI ) * Z;

    if (T1 > 1e-10) CS *= 1.5*T1*std::sqrt(T1);
    return CS;
}

std::vector<double>& ElectronStates::continuous(const std::vector<double>& energy) {

    auto N = new std::vector<double>(energy.size());

    for (std::size_t i = 0; i < energy.size(); ++i) {
        (*N)[i] = continuous(energy[i]);
    }

    return *N;
}

double* ElectronStates::continuous(const double* energy, const std::size_t& n) {

    double* N = new double[n];

    for (std::size_t i = 0; i < n; ++i) {
        N[i] = continuous(energy[i]);
    }

    return N;
}

double ElectronStates::eBoundary() {

    std::vector<double> elvl(nMax*(nMax + 1)/2);
    for (int n = 1; n <= nMax; ++n) {
        auto en = e[n];
        for (int l = 0; l < n; ++l)
            elvl[l + n*(n - 1)/2] = en[l];
    }

    std::sort(elvl.begin(), elvl.end());

    int i = elvl.size() - 2;
    std::vector<double> roots; roots.resize(0);
    while(roots.size() == 0 && i > 1) {
        double dEleft  = 0.5*(elvl[i - 1] - elvl[i - 2]);
        double dEright = 0.5*(elvl[i + 1] - elvl[i]);
        roots = BEroots(elvl[i - 1] - dEleft, elvl[i] + dEright); --i;
    }
    if (roots.size() == 0) return 0.0;
    else return roots[roots.size() - 1];
}

std::vector<double> ElectronStates::BEroots(const double& eLeft, const double& eRight) {
    std::vector<double> roots;
    std::vector<double> eLefts;
    std::vector<double> eRights;
    // search for intervals where sign changes
    int Nintervals = 10;
    std::vector<bool> signs(Nintervals);
    double de = (eRight - eLeft)/(Nintervals - 1);
    for (int i = 0; i < Nintervals; ++i) {
        double e = eLeft + i*de;
        signs[i] = pseudoDS(e) >= pseudoCS(e);
    }
    int nroots = 0;
    for (int i = 0; i < Nintervals - 1; ++i) {
        if (signs[i] != signs[i + 1]) {
            eLefts.push_back(eLeft + i*de);
            eRights.push_back(eLeft + (i + 1)*de);
            ++nroots;
        }
    }
    // calculate roots
    roots.reserve(nroots);
    roots.resize(0);
    for (int i = 0; i < nroots; ++i) {
        double eL = eLefts[i];  double dStatesL = pseudoDS(eL) - pseudoCS(eL);
        double eR = eRights[i]; double dStatesR = pseudoDS(eR) - pseudoCS(eR);
        double eB = 0.5*(eL + eR);
        double error = 1.0;
        while (error > 0.1*tolerance) {
            double dStates = pseudoDS(eB) - pseudoCS(eB);
            if (dStates*dStatesL > 0.0) { eL = eB; dStatesL = dStates; }
            if (dStates*dStatesR > 0.0) { eR = eB; dStatesR = dStates; }
            error = std::abs((eR - eL)/(eR + eL));
            eB = 0.5*(eL + eR);
        }
        if (std::abs(eB - eLefts[i]) > 100*tolerance && 
        	std::abs(eB - eRights[i]) > 100*tolerance) {
            roots.push_back(eB);
        }
    }
    return roots;
}

double ElectronStates::pseudoDS(const double& energy) {
    double DS = 0.0;
    for (int n = 1; n <= nMax; ++n) {
        auto en = e[n];
        for (int l = 0; l < n; ++l) {
            if (en[l] < energy) DS += 2 * l + 1;
        }
    }
    return DS;
}

double ElectronStates::pseudoCS(const double& energy) {
    Array<RHSBE::dim> y;
    y.fill(0.0);
    double Z = V1/VZ;

    RHSBE rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    rhs.set_e(energy*std::pow(Z, -4.0/3.0));

    Solver<PD853<RHSBE>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    solver.integrate(rhs, y, 1.0, 0.0);

    return y[RHSBE::result]*2.0*V1*std::sqrt(2.0) / ( M_PI * M_PI ) * Z;
}
