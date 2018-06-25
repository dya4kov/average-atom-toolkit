#include <cmath>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/ODE/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/ODE/NDV.h>
#include <average-atom-toolkit/thomas-fermi/atom/ODE/NDT.h>
#include <average-atom-toolkit/thomas-fermi/atom/ODE/NDM.h>

#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

#ifdef ENABLE_MULTITHREADING
#include <thread>
#endif

using namespace aatk::TF;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;
using aatk::TF::ODE::RHSPotential;

using aatk::TF::ODE::RHSNDV;
using aatk::TF::ODE::RHSNDT;
using aatk::TF::ODE::RHSNDM;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;

// init static variables
const int    ChemicalPotential::vSize =   181;
const int    ChemicalPotential::tSize =   201;
const double ChemicalPotential::lgV0  = -10.0;
const double ChemicalPotential::lgT0  = -10.0;
const double ChemicalPotential::lgVstep = 0.1;
const double ChemicalPotential::lgTstep = 0.1;
const double ChemicalPotential::bestTolerance  = 1e-12;

ChemicalPotential::ChemicalPotential() : 
#ifdef ENABLE_MULTITHREADING
    threadsLimit(16),
#endif
    tolerance(1e-6), Z(1.0)
{}

ChemicalPotential::ChemicalPotential(const ChemicalPotential& mu) {
    tolerance = mu.tolerance;
    Z = mu.Z;
#ifdef ENABLE_MULTITHREADING
    threadsLimit = mu.threadsLimit;
#endif
}

ChemicalPotential& ChemicalPotential::operator=(const ChemicalPotential& mu) {
    tolerance = mu.tolerance;
    Z = mu.Z;
#ifdef ENABLE_MULTITHREADING
    threadsLimit = mu.threadsLimit;
#endif
    return *this;
}

void ChemicalPotential::setZ(const double& _Z) { Z = _Z; }

#ifdef ENABLE_MULTITHREADING
void ChemicalPotential::setThreadsLimit(const std::size_t& Nthreads) {
    threadsLimit = Nthreads; 
}
#endif

void ChemicalPotential::setTolerance(const double& t) {
    tolerance = std::max(t, bestTolerance);
}

double ChemicalPotential::operator()(const double& V, const double& T) {
    double result;
    bool finished;
    M(V, T, result, finished);
    return result;
}

double ChemicalPotential::DV(const double& V, const double& T) {
    double result;
    bool finished;
    MDV(V, T, result, finished);
    return result;
}

double ChemicalPotential::DT(const double& V, const double& T) {
    double result;
    bool finished;
    MDT(V, T, result, finished);
    return result;
}

double* ChemicalPotential::operator()(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&ChemicalPotential::M, this, _1, _2, _3, _4);
    auto muBuf = evaluate(func, V, T, vsize, tsize);
    double* result = new double[vsize*tsize];
    for (std::size_t i = 0; i < vsize*tsize; ++i) result[i] = muBuf[i]; // copy
    return result;
}

double* ChemicalPotential::DV(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&ChemicalPotential::MDV, this, _1, _2, _3, _4);
    auto muBuf = evaluate(func, V, T, vsize, tsize);
    double* result = new double[vsize*tsize];
    for (std::size_t i = 0; i < vsize*tsize; ++i) result[i] = muBuf[i]; // copy
    return result;
}

double* ChemicalPotential::DT(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&ChemicalPotential::MDT, this, _1, _2, _3, _4);
    auto muBuf = evaluate(func, V, T, vsize, tsize);
    double* result = new double[vsize*tsize];
    for (std::size_t i = 0; i < vsize*tsize; ++i) result[i] = muBuf[i]; // copy
    return result;
}

std::vector<double> ChemicalPotential::operator()(
    const std::vector<double>& V, 
    const std::vector<double>& T 
) {
    auto func = std::bind(&ChemicalPotential::M, this, _1, _2, _3, _4);
    auto result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    return result;
}

std::vector<double> ChemicalPotential::DV(
    const std::vector<double>& V, 
    const std::vector<double>& T 
) {
    auto func = std::bind(&ChemicalPotential::MDV, this, _1, _2, _3, _4);
    auto result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    return result;
}

std::vector<double> ChemicalPotential::DT(
    const std::vector<double>& V, 
    const std::vector<double>& T 
) {
    auto func = std::bind(&ChemicalPotential::MDT, this, _1, _2, _3, _4);
    auto result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    return result;
}

#ifdef ENABLE_MULTITHREADING
void ChemicalPotential::updateThreads(
    std::size_t& threads, 
    std::size_t& current, 
    std::size_t& last, 
    bool* finished
) {
    for (std::size_t thread = current; thread < last; ++thread) {
        if (finished[thread]) {
            --threads; if (threads == 0) break;
        }
    }
    while (finished[current] && current < last) ++current;
}
#endif

std::vector<double> ChemicalPotential::evaluate(
    std::function<void(const double&, const double&, double&, bool&)> func, 
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    std::vector<double> result(vsize*tsize);
    double* c_result = result.data();
    bool*   finished = new bool[vsize*tsize];

#ifdef ENABLE_MULTITHREADING
    std::size_t threads = 0;
    std::size_t current = 0;
    std::size_t last    = 0;
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            finished[v*tsize + t] = false;
            std::thread run(func, std::cref(V[v]), std::cref(T[t]), 
                                  std::ref(c_result[v*tsize + t]), 
                                  std::ref(finished[v*tsize + t]));
            run.detach(); ++threads; ++last;
            while (threads == threadsLimit) {
                updateThreads(threads, current, last, finished);
            }
        }
    }
    bool all_finished = false;
    while (!all_finished) {
        updateThreads(threads, current, last, finished);
        all_finished = (current == last);
    }
#else
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            func(V[v], T[t], c_result[v*tsize + t], finished[v*tsize + t]);
        }
    }
#endif

    delete[] finished;
    return result;
}

void ChemicalPotential::M(
    const double& V, 
    const double& T, 
    double&  result, 
    bool&  finished
) {
    double V1 = V*Z;
    double T1 = T*std::pow(Z, -4.0/3.0);
    double M1 = mu1(V1, T1, tolerance);
    result = M1*std::pow(Z, 4.0/3.0);
    finished = true;
}

void ChemicalPotential::MDV(
    const double& V, 
    const double& T, 
    double&  result, 
    bool&  finished
) {
    double V1 = V*Z;
    double T1 = T*std::pow(Z, -4.0/3.0);
    double M1 = mu1(V1, T1, tolerance);

    Array<RHSNDV::dim> yNDV;
    Array<RHSNDM::dim> yNDM;
    yNDV.fill(0.0);
    yNDM.fill(0.0);

    RHSNDV rhsNDV;     RHSNDM rhsNDM;
    rhsNDV.set_V(V1);  rhsNDM.set_V(V1);
    rhsNDV.set_T(T1);  rhsNDM.set_T(T1);
    rhsNDV.set_mu(M1); rhsNDM.set_mu(M1);

    Solver<PD853<RHSNDV>> solverNDV;
    Solver<PD853<RHSNDM>> solverNDM;
    solverNDV.setTolerance(0.0, 0.1*tolerance);
    solverNDM.setTolerance(0.0, 0.1*tolerance);
    solverNDV.integrate(rhsNDV, yNDV, 1.0, 0.0);
    solverNDM.integrate(rhsNDM, yNDM, 1.0, 0.0);

    double r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);

    double NDV = std::pow(Z, 5.0/3.0)*(r0*yNDV[RHSNDV::result] + 1.0/(3.0*V1));
    double NDM = std::pow(Z, -2.0/3.0)*r0*yNDM[RHSNDM::result];

    result = -NDV/NDM;
    finished = true;
}

void ChemicalPotential::MDT(
    const double& V, 
    const double& T, 
    double&  result, 
    bool&  finished
) {
    double V1 = V*Z;
    double T1 = T*std::pow(Z, -4.0/3.0);
    double M1 = mu1(V1, T1, tolerance);

    Array<RHSNDT::dim> yNDT;
    Array<RHSNDM::dim> yNDM;
    yNDT.fill(0.0);
    yNDM.fill(0.0);

    RHSNDT rhsNDT;     RHSNDM rhsNDM;
    rhsNDT.set_V(V1);  rhsNDM.set_V(V1);
    rhsNDT.set_T(T1);  rhsNDM.set_T(T1);
    rhsNDT.set_mu(M1); rhsNDM.set_mu(M1);

    Solver<PD853<RHSNDT>> solverNDT;
    Solver<PD853<RHSNDM>> solverNDM;
    solverNDT.setTolerance(0.0, 0.1*tolerance);
    solverNDM.setTolerance(0.0, 0.1*tolerance);
    solverNDT.integrate(rhsNDT, yNDT, 1.0, 0.0);
    solverNDM.integrate(rhsNDM, yNDM, 1.0, 0.0);

    double r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);

    double NDT = std::pow(Z, -2.0/3.0)*r0*yNDT[RHSNDT::result];
    double NDM = std::pow(Z, -2.0/3.0)*r0*yNDM[RHSNDM::result];

    result = -NDT/NDM;
    finished = true;
}

double ChemicalPotential::mu1(const double& V1, const double& T1, const double& tol) {
    Array<RHSPotential::dim> phi;

    RHSPotential rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    
    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tol);

    double phi_0 = std::pow(4.0*M_PI/3.0/V1, 1.0/3.0);

    double xFrom = 1.0;
    double xTo   = 0.0;
    double delta = 0.1;

    double muStart;
    if (T1 <= 1e-10) muStart = mu1_approx(std::log10(V1));
    else muStart = mu1_approx(std::log10(V1), std::log10(T1));

    double muPrev = muStart;
    rhs.set_mu(muPrev);
    phi.fill(0.0);
    solver.integrate(rhs, phi, xFrom, xTo);
    double phiPrev = phi[0];

    double muCurr = muStart - delta*std::abs(muStart) - tol;
    rhs.set_mu(muCurr);
    phi.fill(0.0);
    solver.integrate(rhs, phi, xFrom, xTo);
    double phiCurr = phi[0];

    double muNext = 0.0;
    double phiNext = 0.0;

    double error = std::abs(muPrev - muCurr)/std::abs(muPrev + muCurr + tol);

    while (error > tol) {
        muNext = muCurr - (phiCurr - phi_0)*(muCurr - muPrev)/(phiCurr - phiPrev);
        phi.fill(0.0);
        rhs.set_mu(muNext);

        solver.integrate(rhs, phi, xFrom, xTo);
        phiNext = phi[0];

        muPrev  = muCurr;
        phiPrev = phiCurr;
        muCurr  = muNext;
        phiCurr = phiNext;

        error = std::abs(muPrev - muCurr)/std::abs(muPrev + muCurr + tol);
    }

    return muNext;
}


double ChemicalPotential::mu1_approx(const double& lgV) {
    const double B0 =  0.648742997083556;
    const double B1 = -0.704984628856768;
    const double B2 = -0.0224226496439102;
    const double B3 = -0.00419385235723519;
    const double B4 = -3.75915351702641E-5;
    const double B5 =  3.94764845762704E-5;
    const double B6 =  5.4828018180471E-7;
    const double B7 = -1.49964096611993E-7;
    double lgV1 = lgV;
    double lgV2 = lgV1*lgV1;
    double lgV3 = lgV1*lgV2;
    double lgV4 = lgV2*lgV2;
    double lgV5 = lgV3*lgV2;
    double lgV6 = lgV3*lgV3;
    double lgV7 = lgV4*lgV3;
    double lgMu = B0 + B1*lgV1 + B2*lgV2 + B3*lgV3 + B4*lgV4 + B5*lgV5 + B6*lgV6 + B7*lgV7;
    return std::pow(10.0, lgMu);
}

double ChemicalPotential::mu1_approx(const double& lgV, const double& lgT) {
    int v, t;
    v = (int) std::floor((lgV - lgV0)/lgVstep);
    t = (int) std::floor((lgT - lgT0)/lgTstep);
    double result = 0.0;
    if ((v >= 0 && v < vSize) || (t < tSize)) {
        double f00, f01, f10, f11;
        f00 = table[v     +       t*vSize];
        f01 = table[v     + (t + 1)*vSize];
        f10 = table[v + 1 +       t*vSize];
        f11 = table[v + 1 + (t + 1)*vSize];
        int sign = 0;
        sign = (sign == 0 && f00 > 0.0 && f01 > 0.0 && f10 > 0.0 && f11 > 0.0) ?  1 : 0;
        sign = (sign == 0 && f00 < 0.0 && f01 < 0.0 && f10 < 0.0 && f11 < 0.0) ? -1 : 0;
        if (sign != 0) {
            f00 = std::log10(std::abs(f00));
            f01 = std::log10(std::abs(f01));
            f10 = std::log10(std::abs(f10));
            f11 = std::log10(std::abs(f11));
        }
        double V0 = (lgV - lgV0)/lgVstep - 1.0*v;
        double T0 = (lgT - lgT0)/lgTstep - 1.0*t;
        double a[2][2];
        a[0][0] = f00;
        a[1][0] = f10 - f00;
        a[0][1] = f01 - f00;
        a[1][1] = f11 + f00 - (f10 + f01);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                result += a[i][j]*std::pow(V0, i)*std::pow(T0, j);
            }
        }
        if (sign != 0) result = sign*std::pow(10.0, result);
    }
    return result;
}