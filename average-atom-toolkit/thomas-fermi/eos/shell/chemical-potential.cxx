#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/shell/ODE/potential.h>

#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/electron-density.h>

#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/shell/chemical-potential.h>

#ifdef ENABLE_MULTITHREADING
#include <thread>
#endif

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using aatk::TF::shell::ODE::RHSPotential;

using namespace aatk::TF::shell;

#define DRHOSH_INTERPOLATION_SIZE 3501

ChemicalPotential::ChemicalPotential() : 
#ifdef ENABLE_MULTITHREADING
    threadsLimit(std::max(4u, std::thread::hardware_concurrency())), 
#endif
    tolerance(1e-6), Z(1.0), nmax(15)
{}

ChemicalPotential::ChemicalPotential(const ChemicalPotential& dmu) {
    tolerance = dmu.tolerance;
    Z = dmu.Z;
#ifdef ENABLE_MULTITHREADING
    threadsLimit = dmu.threadsLimit;
#endif
}

ChemicalPotential& ChemicalPotential::operator=(const ChemicalPotential& dmu) {
    tolerance = dmu.tolerance;
    Z = dmu.Z;
#ifdef ENABLE_MULTITHREADING
    threadsLimit = dmu.threadsLimit;
#endif
    return *this;
}

void ChemicalPotential::setZ(const double _Z) { Z = _Z; }

#ifdef ENABLE_MULTITHREADING
void ChemicalPotential::setThreadsLimit(const std::size_t Nthreads) {
    threadsLimit = Nthreads; 
}
#endif

void ChemicalPotential::setTolerance(const double t) { tolerance = t; }
void ChemicalPotential::setNmax(const std::size_t _nmax) { nmax = _nmax; }

double ChemicalPotential::operator()(const double V, const double T) { return M(V, T); }

double* ChemicalPotential::operator()(
    const double* V, 
    const double* T, 
    const std::size_t vsize, 
    const std::size_t tsize
) {
    double* result = new double[vsize*tsize];
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            result[v*tsize + t] = M(V[v], T[t]);
        }
    }
    return result;
}

std::vector<double> ChemicalPotential::operator()(
    const std::vector<double>& V, 
    const std::vector<double>& T 
) {
    std::vector<double> result(V.size()*T.size());
    for (std::size_t v = 0; v < V.size(); ++v) {
        for (std::size_t t = 0; t < T.size(); ++t) {
            result[v*T.size() + t] = M(V[v], T[t]);
        }
    }
    return result;
}

double ChemicalPotential::M(
    const double V, 
    const double T
) {

    ::aatk::TF::ElectronStates N;
    N.setVTZ(V, T, Z);
    N.setTolerance(tolerance);
    N.setNmax(nmax);

#ifdef ENABLE_MULTITHREADING
    N.setThreadsLimit(threadsLimit);
#endif
    
    double eb  = N.eBoundary();

    ::aatk::TF::ChemicalPotential M;
    M.setTolerance(tolerance);
    M.setZ(Z);

    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double mu1 = M(V, T)*std::pow(Z, -4.0/3.0);

    // prepare electron density
    std::vector<double> x      (DRHOSH_INTERPOLATION_SIZE, 0.0);
    std::vector<double> sqrtx  (DRHOSH_INTERPOLATION_SIZE, 0.0);

    double dx = 1.0/(x.size() - 1);
    for (std::size_t i = 0; i < x.size(); ++i) {
    	sqrtx[i] = i*dx;
    	x[i] = sqrtx[i]*sqrtx[i];
    }

    ::aatk::TF::shell::ElectronDensity dn;
    dn.setVTZ(V, T, Z);
    dn.setTolerance(tolerance);
    dn.setNmax(nmax);
    dn.setEnergyLevels(N.eLevel());
    dn.setBoundary(eb);

#ifdef ENABLE_MULTITHREADING
    dn.setThreadsLimit(threadsLimit);
#endif

    double r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);

    std::vector<double> dnsh = dn(x);
    dnsh[0] = 0.0;
    for (std::size_t i = 1; i < dnsh.size(); ++i) {
    	dnsh[i] = 4.0*M_PI*r0*r0*x[i]*x[i]*dnsh[i];
    }

    // calculate dmush
    RHSPotential rhs;

    rhs.set_V    (V1);
    rhs.set_T    (T1);
    rhs.set_Z    (Z);
    rhs.set_mu   (mu1);
    rhs.set_eb   (eb);
    rhs.set_drho (sqrtx, dnsh);

    Array<RHSPotential::dim> psi;
    double xFrom = 1.0;
    double xTo   = 0.0;

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(tolerance, 0.0);

    double psiStart  = 0.0;
    double psiPrev_1 = psiStart;
    psi.fill(0.0);
    psi[2] = psi[3]  = psiPrev_1;
    solver.integrate(rhs, psi, xFrom, xTo);
    double psiPrev_0 = psi[2];

    double psiCurr_1 = 0.01*mu1*std::pow(Z, 4.0/3.0);
    psi.fill(0.0);
    psi[2] = psi[3] = psiCurr_1;
    solver.integrate(rhs, psi, xFrom, xTo);
    double psiCurr_0 = psi[2];

    double psiNext_1 = 0.0;
    double psiNext_0 = 0.0;

    double error = std::abs(psiCurr_1 - psiPrev_1)/std::abs(psiCurr_1 + psiPrev_1 + tolerance);
    while (error > tolerance) {
        psiNext_1 = psiCurr_1 - psiCurr_0*(psiCurr_1 - psiPrev_1)/(psiCurr_0 - psiPrev_0);

        psi.fill(0.0);
        psi[2] = psi[3] = psiNext_1;

        solver.integrate(rhs, psi, xFrom, xTo);
        psiNext_0 = psi[2];

        psiPrev_1 = psiCurr_1;
        psiPrev_0 = psiCurr_0;
        psiCurr_1 = psiNext_1;
        psiCurr_0 = psiNext_0;

        // std::cout << "psiCurr_1 = " << psiCurr_1 << ", psiCurr_0 = " << psiCurr_0 << std::endl;

        error = std::abs(psiPrev_1 - psiCurr_1)/std::abs(psiPrev_1 + psiCurr_1 + tolerance);
        // error = std::abs(psiCurr_0);
        // error = std::abs(psiCurr_0);
        // double error2 = std::abs(psiCurr_0);
        // error  = std::sqrt(error1*error1 + error2*error2);
    }

    return psiNext_1;
}