#include <cmath>
#include <thread>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/ODE/potential.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/chemical-potential.h>

using namespace aatk::TF;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;
using aatk::TF::ODE::RHSPotential;

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

// double ChemicalPotential::globalTolerance = 1e-6;
// double ChemicalPotential::globalZ = 1.0;

// std::mutex           ChemicalPotential::mutexBuf;
// std::vector<double>  ChemicalPotential::vBuf;
// std::vector<double>  ChemicalPotential::tBuf;
// std::vector<double>  ChemicalPotential::muBuf;

ChemicalPotential::ChemicalPotential() : 
    tolerance(1e-6), Z(1.0), threadsLimit(4)
{}

ChemicalPotential::ChemicalPotential(const ChemicalPotential& mu) {
    tolerance = mu.tolerance;
    Z = mu.Z;
    threadsLimit = mu.threadsLimit;
}

ChemicalPotential& ChemicalPotential::operator=(const ChemicalPotential& mu) {
    tolerance = mu.tolerance;
    Z = mu.Z;
    threadsLimit = mu.threadsLimit;
    return *this;
}

void ChemicalPotential::setZ(const double& _Z) { Z = _Z; }

void ChemicalPotential::setThreadsLimit(const std::size_t& Nthreads) {
    threadsLimit = Nthreads; 
}

void ChemicalPotential::setTolerance(const double& t) {
    tolerance = std::max(t, bestTolerance);
}

double ChemicalPotential::operator()(const double& V, const double& T) {
    double result;
    double V1 = V*Z;
    double T1 = T*std::pow(Z, -4.0/3.0);
    result = mu1(V1, T1, tolerance);
    return result*std::pow(Z, 4.0/3.0);
}

double* ChemicalPotential::operator()(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&ChemicalPotential::mu, this, _1, _2, _3, _4);
    auto muBuf = evaluate(func, V, T, vsize, tsize);
    double* result = new double[vsize*tsize];
    for (std::size_t i = 0; i < vsize*tsize; ++i) result[i] = muBuf[i]; // copy
    return result;
}

std::vector<double> ChemicalPotential::operator()(
    const std::vector<double>& V, 
    const std::vector<double>& T 
) {
    auto func = std::bind(&ChemicalPotential::mu, this, _1, _2, _3, _4);
    auto result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    return result;
}

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
                for (std::size_t thread = current; thread < last; ++thread) {
                    if (finished[thread]) --threads;
                }
                while (finished[current] && current < last) ++current;
            }
        }
    }
    bool all_finished = false;
    while (!all_finished) {
        while (finished[current] && current < last) ++current;
        all_finished = (current == last);
    }
    delete[] finished;
    return result;
}

void ChemicalPotential::mu(
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

// bool ChemicalPotential::validBuf() {
//     bool valid = true;
//     valid = valid && vBuf.size() > 1;
//     valid = valid && tBuf.size() > 1;
//     valid = valid && (std::abs(localTolerance - globalTolerance) < bestTolerance);
//     valid = valid && (std::abs(localZ - globalZ) < bestTolerance);
//     return valid;
// }

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

    bool convergenceSuccess = false;
    double mu;
    while (!convergenceSuccess) {
        double muLeftStart  = muStart - delta*std::max(10.0, std::abs(muStart));
        double muRightStart = muStart + delta*std::max(10.0, std::abs(muStart));
        double muLeft  = muLeftStart;
        double muRight = muRightStart;
        double error   = std::abs(muRight - muLeft)/std::abs(muLeft + muRight + tol);
        while (error > tol) {
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
        mu = 0.5*(muLeft + muRight);
        convergenceSuccess = true;
        if (std::abs( muLeftStart - muLeft )/std::abs( muLeftStart + muLeft ) < tol*100.0 || 
            std::abs(muRightStart - muRight)/std::abs(muRightStart + muRight) < tol*100.0  ) {
            convergenceSuccess = false;
            delta = delta + 0.1;
            if (delta > 0.45) {
                std::cout << "ChemicalPotential error: too big delta" << std::endl;
                exit(0);
            }
        }
    }
    return mu;
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

// int ChemicalPotential::vIndex(const double& V) {
//     int size = vBuf.size();
//     int ileft  = 0;
//     int iright = size;
//     int i = (ileft + iright)/2;
//     while (std::abs(vBuf[i] - V) > bestTolerance) {
//         if (V > vBuf[i]) ileft = i;
//         else iright = i;
//         if (ileft == iright) break;
//         i = (ileft + iright)/2;
//     }
//     return ileft == iright ? -1 : i;
// }
// 
// int ChemicalPotential::tIndex(const double& T) {
//     int size = tBuf.size();
//     int ileft  = 0;
//     int iright = size;
//     int i = (ileft + iright)/2;
//     while (std::abs(tBuf[i] - T) > bestTolerance) {
//         if (T > tBuf[i]) ileft = i;
//         else iright = i;
//         if (ileft == iright) break;
//         i = (ileft + iright)/2;
//     }
//     return ileft == iright ? -1 : i;
// }
// 
// void ChemicalPotential::prepareBuf(
//     const double* V, 
//     const double* T, 
//     const std::size_t& vsize, 
//     const std::size_t& tsize
// ) {
//     std::lock_guard<std::mutex> lock(mutexBuf);
//     globalZ = localZ;
//     globalTolerance = localTolerance;
//     vBuf = sorted(V, vsize);
//     tBuf = sorted(T, tsize);
//     muBuf.resize(vBuf.size()*tBuf.size());
// }
// 
// std::vector<double> ChemicalPotential::sorted(
//     const double*     array,
//     const std::size_t& size
// ) {
//     std::vector<double> result(array, array + size);
//     std::sort(result.begin(), result.end());
//     return result;
// }
