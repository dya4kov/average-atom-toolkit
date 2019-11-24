#include <cmath>
#include <numeric>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/wave-function.h>
#include <average-atom-toolkit/thomas-fermi/atom/electron-density.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/electron-density.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

#ifdef ENABLE_MULTITHREADING
#include <thread>
#endif

using namespace aatk::TF::shell;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::std::placeholders::_4;

ElectronDensity::ElectronDensity() :
    V(1.0), T(1.0), Z(1.0), mu(4.100577730112),
    tolerance(1e-6), eBoundary(1e+15), 
    nmax(15), 
#ifdef ENABLE_MULTITHREADING
    threadsLimit(std::max(4u, std::thread::hardware_concurrency())),
#else
    threadsLimit(0),
#endif
    ready(false)
{
    p_accumulate = std::bind(&ElectronDensity::accumulate, this, _1, _2, _3, _4);
}

ElectronDensity::ElectronDensity(const ElectronDensity& eDens) {
    V            = eDens.V; 
    T            = eDens.T; 
    Z            = eDens.Z; 
    mu           = eDens.mu;
    tolerance    = eDens.tolerance;
    eBoundary    = eDens.eBoundary;
    e            = eDens.e;
    nmax         = eDens.nmax;
    ready        = eDens.ready;
    threadsLimit = eDens.threadsLimit;
    p_accumulate = std::bind(&ElectronDensity::accumulate, this, _1, _2, _3, _4);
}

ElectronDensity& ElectronDensity::operator=(const ElectronDensity& eDens) {
    V            = eDens.V; 
    T            = eDens.T; 
    Z            = eDens.Z; 
    mu           = eDens.mu;
    tolerance    = eDens.tolerance;
    eBoundary    = eDens.eBoundary;
    e            = eDens.e;
    nmax         = eDens.nmax;
    ready        = eDens.ready;
    threadsLimit = eDens.threadsLimit;
    p_accumulate = std::bind(&ElectronDensity::accumulate, this, _1, _2, _3, _4);
    return *this;
}

void ElectronDensity::setEnergyLevels(const ::aatk::TF::EnergyLevel& _e) {
    e = _e;
}

void ElectronDensity::setTolerance(const double t) {
    tolerance = t;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setTolerance(t);
    ready = false;
}

#ifdef ENABLE_MULTITHREADING
void ElectronDensity::setThreadsLimit(const std::size_t Nthreads) {
    threadsLimit = Nthreads;
    e.setThreadsLimit(Nthreads);
}
#endif

void ElectronDensity::setNmax(const std::size_t Nmax) { nmax = Nmax; }

void ElectronDensity::setBoundary(const double eb) { eBoundary = eb; }

void ElectronDensity::setV(const double _V) {
    V = _V;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setV(V);
    ready = false;
}

void ElectronDensity::setT(const double _T) {
    T = _T;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setT(T);
    ready = false;
}

void ElectronDensity::setZ(const double _Z) {
    Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setZ(Z);
    ready = false;
}

void ElectronDensity::setVTZ(
    const double _V,
    const double _T,
    const double _Z
) {
    V = _V; T = _T; Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(tolerance);
    mu = M(V, T);
    e.setVTZ(V, T, Z);
    ready = false;
}

double ElectronDensity::operator()(const double x) {

    if (!ready) e.prepareLevelsBelow(nmax);
    ready = true;

    double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);

    ::aatk::TF::shell::WaveFunction WF;
    WF.setVTZ(V, T, Z);
    WF.setTolerance(tolerance);

    ::aatk::TF::ElectronDensity rhotf;
    rhotf.setVTZ(V, T, Z);
    rhotf.setTolerance(tolerance);
    rhotf.setBoundary(eBoundary);

    double rho = 0.0;
    if (x > 0.0) {
        for (int n = 1; n < nmax; ++n) {
            for (int l = 0; l < n; ++l) {
                double enl    = e(n,l);
                if (enl > eBoundary) continue;
                double Nnl    = (2.0*l + 1)/(1.0 + std::exp((enl - mu)/T));
                double lambda = l + 0.5;
                double Rnl    = WF(enl, lambda, x);
                rho          += Nnl*Rnl*Rnl;
            }
        }

        r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
        rho *= 1.0/(2.0*M_PI*r0*r0*x*x);
    }

    return rho - rhotf(x);
}

void ElectronDensity::accumulate(const double* x, double* rho, const std::size_t size, const std::size_t ithread) {
    ::aatk::TF::shell::WaveFunction Rnl;
    Rnl.setVTZ(V, T, Z);
    Rnl.setTolerance(tolerance);

    std::vector<double> wf(size);
    std::size_t threadsCounter = 0;

    double T1 = T*std::pow(Z, -4.0/3.0);

    for (std::size_t n = 1; n <= nmax; ++n) {
        for (std::size_t l = 0; l < n; ++l) {
            if (ithread == threadsCounter) {
                double enl = e(n,l);
                if (enl <= eBoundary) {
                    double lambda = l + 0.5;
                    Rnl(enl, lambda, x, wf.data(), size); // set wf
                    double Nnl;
                    if (T1 > 1e-10)
                         Nnl = (2.0*l + 1.0)/(1.0 + std::exp((enl - mu)/T));
                    else Nnl = enl <= mu ? 2.0*l + 1.0 : 0.0;
                    for (std::size_t i = 0; i < size; ++i) {
                        if (x[i] > 0.0)
                            rho[i] += 1.0/(x[i]*x[i])*Nnl*wf[i]*wf[i];
                    }
                }
            }
            if (threadsLimit > 0)
                threadsCounter = (threadsCounter + 1) % threadsLimit;
        }
    }

    return;
}

void ElectronDensity::operator()(const double* x, double* rho, const std::size_t size) {
    if (!ready) e.prepareLevelsBelow(nmax);
    ready = true;

    for (std::size_t i = 0; i < size; ++i) rho[i] = 0.0;

#ifdef ENABLE_MULTITHREADING
    std::vector<std::vector<double>> rho_(threadsLimit, std::vector<double>(size, 0.0));
    
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_accumulate, x, rho_[ithread].data(), size, ithread));
    }
    for (auto&& thread : threads) thread.join();

    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        for (std::size_t i = 0; i < size; ++i) {
            rho[i] += rho_[ithread][i];
        }
    }
#else
    accumulate(x, rho, size, 0);
#endif

    ::aatk::TF::ElectronDensity rhotf;
    rhotf.setVTZ(V, T, Z);
    rhotf.setTolerance(tolerance);
    rhotf.setBoundary(eBoundary);
    auto rhotf_calc = rhotf(x, size);

    double r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
    for (std::size_t i = 0; i < size; ++i) {
        rho[i] *= 1.0/(2.0*M_PI*r0*r0);
        rho[i] -= rhotf_calc[i];
    }
}

double* ElectronDensity::operator()(const double* x, const std::size_t size) {
    double* rho = new double[size];
    operator()(x, rho, size);
    return rho;
}

std::vector<double> ElectronDensity::operator()(const std::vector<double>& x) {
    std::vector<double> rho(x.size());
    operator()(x.data(), rho.data(), x.size());
    return rho;
}