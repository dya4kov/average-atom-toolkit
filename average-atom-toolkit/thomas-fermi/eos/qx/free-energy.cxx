#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/Yfunction.h>

#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/qx/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/qx/free-energy.h>
#include <average-atom-toolkit/thomas-fermi/eos/qx/ODE/F.h>
#include <average-atom-toolkit/thomas-fermi/eos/qx/ODE/FDT.h>

#ifdef ENABLE_MULTITHREADING
#include <thread>
#endif

using numtk::ODE::Array;
using numtk::ODE::Dimension;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using numtk::specfunc::FermiDirac;
using numtk::specfunc::Yfunction;
using numtk::specfunc::FD::ThreeHalf;
using numtk::specfunc::FD::Half;
using numtk::specfunc::FD::MHalf;

using aatk::TF::qx::ODE::RHSdF;
using aatk::TF::qx::ODE::RHSdS;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;
using std::placeholders::_5;
using std::placeholders::_6;
using std::placeholders::_7;

using namespace aatk::TF::qx;

FreeEnergy::FreeEnergy() : 
#ifdef ENABLE_MULTITHREADING
    threadsLimit(std::max(4u, std::thread::hardware_concurrency())),
#endif
    tolerance(1e-6),
    Z(1.0), 
    dE0(0.26990017) 
{
    p_evaluate = std::bind(&FreeEnergy::evaluate, this, _1, _2, _3, _4, _5, _6, _7);
}

double FreeEnergy::operator() (const double V, const double T) { return F    (V, T); }
double FreeEnergy::DV         (const double V, const double T) { return FDV  (V, T); }
double FreeEnergy::DT         (const double V, const double T) { return FDT  (V, T); }
double FreeEnergy::D2V        (const double V, const double T) { return FD2V (V, T); }
double FreeEnergy::DVT        (const double V, const double T) { return FDVT (V, T); }
double FreeEnergy::D2T        (const double V, const double T) { return FD2T (V, T); }

double* FreeEnergy::operator()(
    const double* V,
    const double* T, 
    const std::size_t vsize, 
    const std::size_t tsize
) {
    auto func = std::bind(&FreeEnergy::F, this, _1, _2);
    double* result = new double[vsize*tsize];
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V, T, result, vsize, tsize, ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V, T, result, vsize, tsize, 0);
#endif
    return result;
}

double* FreeEnergy::DV(
    const double* V,
    const double* T, 
    const std::size_t vsize, 
    const std::size_t tsize
) {
    auto func = std::bind(&FreeEnergy::FDV, this, _1, _2);
    double* result = new double[vsize*tsize];
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V, T, result, vsize, tsize, ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V, T, result, vsize, tsize, 0);
#endif
    return result;
}

double* FreeEnergy::DT(
    const double* V,
    const double* T, 
    const std::size_t vsize, 
    const std::size_t tsize
) {
    auto func = std::bind(&FreeEnergy::FDT, this, _1, _2);
    double* result = new double[vsize*tsize];
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V, T, result, vsize, tsize, ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V, T, result, vsize, tsize, 0);
#endif
    return result;
}

double* FreeEnergy::D2V(
    const double* V, 
    const double* T, 
    const std::size_t vsize, 
    const std::size_t tsize
) {
    auto func = std::bind(&FreeEnergy::FD2V, this, _1, _2);
    double* result = new double[vsize*tsize];
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V, T, result, vsize, tsize, ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V, T, result, vsize, tsize, 0);
#endif
    return result;
}

double* FreeEnergy::DVT(
    const double* V, 
    const double* T, 
    const std::size_t vsize, 
    const std::size_t tsize
) {
    auto func = std::bind(&FreeEnergy::FDVT, this, _1, _2);
    double* result = new double[vsize*tsize];
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V, T, result, vsize, tsize, ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V, T, result, vsize, tsize, 0);
#endif
    return result;
}

double* FreeEnergy::D2T(
    const double* V,
    const double* T, 
    const std::size_t vsize, 
    const std::size_t tsize
) {
    auto func = std::bind(&FreeEnergy::FD2T, this, _1, _2);
    double* result = new double[vsize*tsize];
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V, T, result, vsize, tsize, ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V, T, result, vsize, tsize, 0);
#endif
    return result;
}

std::vector<double> FreeEnergy::operator() (
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func = std::bind(&FreeEnergy::F, this, _1, _2);
    std::vector<double> result(V.size()*T.size());
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V.data(), T.data(), result.data(), V.size(), T.size(), ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V.data(), T.data(), result.data(), V.size(), T.size(), 0);
#endif
    return result;
}

std::vector<double> FreeEnergy::DV(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func = std::bind(&FreeEnergy::FDV, this, _1, _2);
    std::vector<double> result(V.size()*T.size());
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V.data(), T.data(), result.data(), V.size(), T.size(), ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V.data(), T.data(), result.data(), V.size(), T.size(), 0);
#endif
    return result;
}

std::vector<double> FreeEnergy::DT(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func = std::bind(&FreeEnergy::FDT, this, _1, _2);
    std::vector<double> result(V.size()*T.size());
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V.data(), T.data(), result.data(), V.size(), T.size(), ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V.data(), T.data(), result.data(), V.size(), T.size(), 0);
#endif
    return result;
}

std::vector<double> FreeEnergy::D2V(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func = std::bind(&FreeEnergy::FD2V, this, _1, _2);
    std::vector<double> result(V.size()*T.size());
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V.data(), T.data(), result.data(), V.size(), T.size(), ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V.data(), T.data(), result.data(), V.size(), T.size(), 0);
#endif
    return result;
}

std::vector<double> FreeEnergy::DVT(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func = std::bind(&FreeEnergy::FDVT, this, _1, _2);
    std::vector<double> result(V.size()*T.size());
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V.data(), T.data(), result.data(), V.size(), T.size(), ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V.data(), T.data(), result.data(), V.size(), T.size(), 0);
#endif
    return result;
}

std::vector<double> FreeEnergy::D2T(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func = std::bind(&FreeEnergy::FD2T, this, _1, _2);
    std::vector<double> result(V.size()*T.size());
#ifdef ENABLE_MULTITHREADING
    std::vector<std::thread> threads;
    for (std::size_t ithread = 0; ithread < threadsLimit; ++ithread) {
        threads.push_back(std::thread(p_evaluate, func, V.data(), T.data(), result.data(), V.size(), T.size(), ithread));
    }
    for (auto&& thread : threads) thread.join();
#else
    evaluate(func, V.data(), T.data(), result.data(), V.size(), T.size(), 0);
#endif
    return result;
}

void FreeEnergy::setTolerance(const double eps) { tolerance = eps; }
void FreeEnergy::setZ(const double _Z) { Z = _Z; }

#ifdef ENABLE_MULTITHREADING
void FreeEnergy::setThreadsLimit(const std::size_t N) {
    threadsLimit = std::max(1LU, N);
}
#endif

void FreeEnergy::evaluate(
    std::function<double(const double, const double)> func, 
    const double* V, 
    const double* T, 
          double* result,
    const std::size_t vsize, 
    const std::size_t tsize,
    const std::size_t ithread
) {
#ifdef ENABLE_MULTITHREADING
    std::size_t max_index = vsize*tsize;
    std::size_t index = ithread;
    while (index < max_index) {
        std::size_t v = index / tsize;
        std::size_t t = index % tsize;
        result[index] = func(V[v], T[t]);
        index += threadsLimit;
    }
#else
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            result[v*tsize + t] = func(V[v], T[t]);
        }
    }
#endif
}

double FreeEnergy::F(const double V, const double T) {
    
    ::aatk::TF::ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    ::aatk::TF::qx::ChemicalPotential dmu;
    dmu.setZ(Z);
    dmu.setTolerance(tolerance);

    double   M1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  dM1 = dmu(V, T)*std::pow(Z, -2.0/3.0);
    double   V1 = V*Z;
    double   T1 = T*std::pow(Z, -4.0/3.0);
    double psi1;

    FermiDirac<MHalf> FDmhalf;

    if (T1 <= 1e-10) psi1 = (6.0 * M_PI)/std::sqrt(2.0)*dM1 - std::sqrt(M1);
    else psi1 = (6.0 * M_PI)/std::sqrt(2.0)*dM1 - 0.5*std::sqrt(T1)*FDmhalf(M1/T1);

    RHSdF rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(M1);

    Array<RHSdF::dim> dF;

    double xFrom = 1.0;
    double xTo   = 0.0;

    dF[0] = dF[1] = 0.0;
    dF[2] = dF[3] = psi1;

    Solver<PD853<RHSdF>> solver;

    solver.setTolerance(0.0, 0.1*tolerance);

    solver.integrate(rhs, dF, xFrom, xTo);

    dF[RHSdF::result] *= rhs.param();
    dF[RHSdF::result] += dE0;
    dF[RHSdF::result] *= std::pow(Z, 5.0/3.0);

    return dF[RHSdF::result];
}

double FreeEnergy::FDV(const double V, const double T) {

    ::aatk::TF::ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    ::aatk::TF::qx::ChemicalPotential dmu;
    dmu.setZ(Z);
    dmu.setTolerance(tolerance);

    double   M1 =  mu(V, T)*std::pow(Z, -4.0/3.0);
    double  dM1 = dmu(V, T)*std::pow(Z, -2.0/3.0);
    double   V1 = V*Z;
    double   T1 = T*std::pow(Z, -4.0/3.0);
    double psi1;

    FermiDirac<MHalf> FDmhalf;

    if (T1 <= 1e-10) psi1 = (6.0 * M_PI)/std::sqrt(2.0)*dM1 - std::sqrt(M1);
    else psi1 = (6.0 * M_PI)/std::sqrt(2.0)*dM1 - 0.5*std::sqrt(T1)*FDmhalf(M1/T1);

    FermiDirac<Half> FDhalf;
    Yfunction Y;

    double result;
    if (T > 1e-10) {
        result = -T1*std::sqrt(T1)/(3.0*M_PI*M_PI*M_PI)*(FDhalf(M1/T1)*psi1 + std::sqrt(T1)*Y(M1/T1));
    }
    else {
        result = -std::sqrt(M1)*M1/(9.0*M_PI*M_PI*M_PI)*(2.0*psi1 + 11.0*std::sqrt(M1));
    }

    result *= std::pow(Z, 8.0/3.0); 
    return result;
}

double FreeEnergy::FDT(const double V, const double T) {

    ::aatk::TF::ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    ::aatk::TF::qx::ChemicalPotential dmu;
    dmu.setZ(Z);
    dmu.setTolerance(tolerance);

    double   M1 =  mu(V, T)*std::pow(Z, -4.0/3.0);
    double  dM1 = dmu(V, T)*std::pow(Z, -2.0/3.0);
    double   V1 = V*Z;
    double   T1 = T*std::pow(Z, -4.0/3.0);
    double psi1;

    FermiDirac<MHalf> FDmhalf;

    if (T1 <= 1e-10) psi1 = (6.0 * M_PI)/std::sqrt(2.0)*dM1 - std::sqrt(M1);
    else psi1 = (6.0 * M_PI)/std::sqrt(2.0)*dM1 - 0.5*std::sqrt(T1)*FDmhalf(M1/T1);

    RHSdS rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(M1);

    Array<RHSdS::dim> dS;

    double xFrom  = 1.0;
    double xTo    = 0.0;
    dS[0] = dS[1] = 0.0;
    dS[2] = dS[3] = psi1;

    Solver<PD853<RHSdS>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    if (T > 1e-10) solver.integrate(rhs, dS, xFrom, xTo);

    double dpsi_0 = dS[3];

    dS[RHSdS::result] *= rhs.param();
    dS[RHSdS::result] += std::sqrt(2.0)/(6.0*M_PI*T1)*dpsi_0;
    dS[RHSdS::result] *= std::pow(Z, 1.0/3.0);

    return -dS[RHSdS::result];
}

double FreeEnergy::FD2V(const double V, const double T) {

    double dV = std::sqrt(std::sqrt(tolerance))*V;

    double FDVleft2  = FDV(V - 2*dV, T);
    double FDVleft1  = FDV(V -   dV, T);
    double FDVright1 = FDV(V +   dV, T);
    double FDVright2 = FDV(V + 2*dV, T);

    return (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dV);
}

double FreeEnergy::FDVT(const double V, const double T) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    double FDVleft2  = FDV(V, T - 2*dT);
    double FDVleft1  = FDV(V, T -   dT);
    double FDVright1 = FDV(V, T +   dT);
    double FDVright2 = FDV(V, T + 2*dT);

    return (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dT);
}

double FreeEnergy::FD2T(const double V, const double T) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    double FDTleft2  = FDT(V, T - 2*dT);
    double FDTleft1  = FDT(V, T -   dT);
    double FDTright1 = FDT(V, T +   dT);
    double FDTright2 = FDT(V, T + 2*dT);

    return (-FDTright2 + 8*FDTright1 - 8*FDTleft1 + FDTleft2)/(12.0*dT);
}