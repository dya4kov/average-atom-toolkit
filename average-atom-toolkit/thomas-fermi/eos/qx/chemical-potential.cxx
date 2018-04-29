#include <cmath>
#include <thread>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

#include <average-atom-toolkit/thomas-fermi/atom/qx/ODE/potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/qx/chemical-potential.h>

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::MHalf;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;

using ::numtk::ODE::Array;
using ::numtk::ODE::Solver;
using ::numtk::ODE::stepper::PD853;

using ::aatk::TF::qx::ODE::RHSPotential;

using namespace aatk::TF::qx;

// init static variables
const int    ChemicalPotential::vSize =   181;
const int    ChemicalPotential::tSize =   201;
const double ChemicalPotential::lgV0  = -10.0;
const double ChemicalPotential::lgT0  = -10.0;
const double ChemicalPotential::lgVstep = 0.1;
const double ChemicalPotential::lgTstep = 0.1;
const double ChemicalPotential::bestTolerance  = 1e-12;

ChemicalPotential::ChemicalPotential() : 
    tolerance(1e-6), Z(1.0), threadsLimit(8)
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
    tolerance = t;
}

double ChemicalPotential::operator()(const double& V, const double& T) {
	double result;
	bool finished;
    M(V, T, result, finished);
    return result;
}

double* ChemicalPotential::operator()(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&ChemicalPotential::M, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

std::vector<double> ChemicalPotential::operator()(
    const std::vector<double>& V, 
    const std::vector<double>& T 
) {
    auto func = std::bind(&ChemicalPotential::M, this, _1, _2, _3, _4);
    double* result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    std::vector<double> vresult(result, result + V.size()*T.size());
    delete[] result;
    return vresult;
}

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

double* ChemicalPotential::evaluate(
    std::function<void(const double&, const double&, double&, bool&)> func, 
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    double* result = new double[vsize*tsize];
    bool* finished = new bool[vsize*tsize];
    std::size_t threads = 0;
    std::size_t current = 0;
    std::size_t last    = 0;
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            finished[v*tsize + t] = false;
            std::thread run(func, std::cref(V[v]), std::cref(T[t]), 
                                  std::ref(result[v*tsize + t]), 
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
    delete[] finished;
    return result;
}

void ChemicalPotential::M(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

    ::aatk::TF::ChemicalPotential mu;
    mu.setZ(Z);    
    mu.setTolerance(tolerance);

    double M1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double V1 = V*Z;
    double T1 = T*std::pow(Z, -4.0/3.0);

    double dM1 = mu1(V1, T1, M1, tolerance);
    result = dM1*std::pow(Z, 2.0/3.0);
    finished = true;
}

double ChemicalPotential::mu1(const double& V1, const double& T1, const double& M1, const double& tol) {

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tol);

    Array<RHSPotential::dim> psi;
    RHSPotential rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(M1);

    double xFrom = 1.0;
    double xTo   = 0.0;
    double delta = 0.1;

    double psiStart;
    if (T1 <= 1e-10) psiStart = mu1_approx(std::log10(V1));
    else psiStart = mu1_approx(std::log10(V1), std::log10(T1));

    double psiPrev_1 = psiStart;
    psi.fill(0.0);
    psi[2] = psi[3] = psiPrev_1;
    solver.integrate(rhs, psi, xFrom, xTo);
    double psiPrev_0 = psi[2];

    double psiCurr_1 = psiStart - delta*std::abs(psiStart) - tol;
    psi.fill(0.0);
    psi[2] = psi[3] = psiCurr_1;
    solver.integrate(rhs, psi, xFrom, xTo);
    double psiCurr_0 = psi[2];

    double psiNext_1 = 0.0;
    double psiNext_0 = 0.0;

    double error = std::abs(psiCurr_1 - psiPrev_1)/std::abs(psiCurr_1 + psiPrev_1 + tol);
    while (error > tol) {
        psiNext_1 = psiCurr_1 - psiCurr_0*(psiCurr_1 - psiPrev_1)/(psiCurr_0 - psiPrev_0);

        psi.fill(0.0);
        psi[2] = psi[3] = psiNext_1;

        solver.integrate(rhs, psi, xFrom, xTo);
        psiNext_0 = psi[2];

        psiPrev_1 = psiCurr_1;
        psiPrev_0 = psiCurr_0;
        psiCurr_1 = psiNext_1;
        psiCurr_0 = psiNext_0;

        error = std::abs(psiPrev_1 - psiCurr_1)/std::abs(psiPrev_1 + psiCurr_1 + tol);
    }

    double psi_1  = psiNext_1;
    double result = std::sqrt(2.0)/ (6.0 * M_PI);
    FermiDirac<MHalf> FDmhalf;

    if (T1 <= 1e-10) result *= std::sqrt(M1) + psi_1;
    else result *= 0.5*std::sqrt(T1)*FDmhalf(M1/T1) + psi_1;

    return result;
}

double ChemicalPotential::mu1_approx(const double& lgV, const double& lgT) {
    int v, t;
    v = (int) std::floor((lgV - lgV0)/lgVstep);
    t = (int) std::floor((lgT - lgT0)/lgTstep);
    double result = 0.0;
    if ((v >= 0 && v < vSize) || (t < tSize)) {
        double f00, f01, f10, f11;
        f00 = -table[v     +       t*vSize];
        f01 = -table[v     + (t + 1)*vSize];
        f10 = -table[v + 1 +       t*vSize];
        f11 = -table[v + 1 + (t + 1)*vSize];
        f00 = std::log10(std::abs(f00));
        f01 = std::log10(std::abs(f01));
        f10 = std::log10(std::abs(f10));
        f11 = std::log10(std::abs(f11));
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
        result = -std::pow(10.0, result);
    }
    return result;
}

double ChemicalPotential::mu1_approx(const double& lgV) {
    const double B0 =  1.20752290577393;
    const double B1 = -0.338439954152781;
    const double B2 = -0.00667391360555085;
    const double B3 = -0.00202025760222042;
    const double B4 = -1.77130362887407E-4;
    const double B5 =  8.67013549505259E-6;
    const double B6 =  1.89781271548473E-6;
    const double B7 =  6.62246701639077E-8;
    double lgV1 = lgV;
    double lgV2 = lgV1*lgV1;
    double lgV3 = lgV1*lgV2;
    double lgV4 = lgV2*lgV2;
    double lgV5 = lgV3*lgV2;
    double lgV6 = lgV3*lgV3;
    double lgV7 = lgV4*lgV3;
    double lgPsi = B0 + B1*lgV1 + B2*lgV2 + B3*lgV3 + B4*lgV4 + B5*lgV5 + B6*lgV6 + B7*lgV7;
    return -std::pow(10.0, lgPsi);
}