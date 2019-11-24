#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/eos/free-energy.h>
#include <average-atom-toolkit/thomas-fermi/eos/ODE/F.h>
#include <average-atom-toolkit/thomas-fermi/eos/ODE/FDT.h>
#include <average-atom-toolkit/thomas-fermi/eos/ODE/FD2T.h>

#ifdef ENABLE_MULTITHREADING
#include <thread>
#endif

using numtk::ODE::Array;
using numtk::ODE::Dimension;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::ThreeHalf;
using ::numtk::specfunc::FD::Half;

using aatk::TF::ODE::RHSF1;
using aatk::TF::ODE::RHSF2;
using aatk::TF::ODE::RHSFDT;
using aatk::TF::ODE::RHSFD2T;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;
using std::placeholders::_5;
using std::placeholders::_6;
using std::placeholders::_7;

using namespace aatk::TF;

FreeEnergy::FreeEnergy() : 
#ifdef ENABLE_MULTITHREADING
  threadsLimit(std::max(4u, std::thread::hardware_concurrency())),
#endif
  tolerance(1e-6), Z(1.0), 
  E0(0.76874512422) 
{
    p_evaluate = std::bind(&FreeEnergy::evaluate, this, _1, _2, _3, _4, _5, _6, _7);
}

double FreeEnergy::operator() (const double V, const double T) { return F    (V, T); }
double FreeEnergy::DV         (const double V, const double T) { return FDV  (V, T); }
double FreeEnergy::DT         (const double V, const double T) { return FDT  (V, T); }
double FreeEnergy::D2V        (const double V, const double T) { return FD2V (V, T); }
double FreeEnergy::D2V2       (const double V, const double T) { return FD2V2(V, T); }
double FreeEnergy::DVT        (const double V, const double T) { return FDVT (V, T); }
double FreeEnergy::DVT2       (const double V, const double T) { return FDVT2(V, T); }
double FreeEnergy::D2T        (const double V, const double T) { return FD2T (V, T); }
double FreeEnergy::D2T2       (const double V, const double T) { return FD2T2(V, T); }

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

    double result;

    ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    RHSF1 rhs1;       RHSF2 rhs2;
    rhs1.set_V(V1);   rhs2.set_V(V1);
    rhs1.set_T(T1);   rhs2.set_T(T1);
    rhs1.set_mu(mu1); rhs2.set_mu(mu1);

    Array<RHSF1::dim> F1;
    Array<RHSF2::dim> F2;

    double xFrom = 1.0;
    double xTo   = 0.0;
    F1.fill(0.0);
    F2.fill(0.0);

    Solver<PD853<RHSF1>> solver1;
    Solver<PD853<RHSF2>> solver2;
    solver1.setTolerance(0.0, 0.1*tolerance);
    solver2.setTolerance(0.0, 0.1*tolerance);

    solver1.integrate(rhs1, F1, xFrom, xTo);
    solver2.integrate(rhs2, F2, xFrom, xTo);

    result = E0 + F1[RHSF1::result]*rhs1.param()
                + F2[RHSF2::result]*rhs2.param();

    result += 0.5*mu1;
    result *= std::pow(Z, 7.0/3.0);

    return result;

}

double FreeEnergy::FDV(const double V, const double T) {

    ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    FermiDirac<ThreeHalf> FD3half;

    double FDV;
    if (T > 1e-10) FDV = -std::pow(2.0*T, 2.5)/(6.0*M_PI*M_PI)*FD3half(mu(V, T)/T);
    else FDV = -std::pow(2.0*mu(V, T), 2.5)/(15.0*M_PI*M_PI);

    return FDV;
}

double FreeEnergy::FDT(const double V, const double T) {

    ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    RHSFDT rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSFDT::dim> FDT;

    double xFrom = 1.0;
    double xTo   = 0.0;
    FDT.fill(0.0);

    Solver<PD853<RHSFDT>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    if (T > 1e-10) solver.integrate(rhs, FDT, xFrom, xTo);

    return Z*FDT[RHSFDT::result]*rhs.param();
}

double FreeEnergy::FD2V(const double V, const double T) {

    double dV = std::sqrt(std::sqrt(tolerance))*V;

    double FDVleft2  = FDV(V - 2*dV, T);
    double FDVleft1  = FDV(V -   dV, T);
    double FDVright1 = FDV(V +   dV, T);
    double FDVright2 = FDV(V + 2*dV, T);

    return (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dV);
}

double FreeEnergy::FD2V2(const double V, const double T) {

    ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    FermiDirac<Half> FDhalf;

    double FD2V;
    if (T > 1e-10) FD2V = -std::sqrt(2.0*T)*T/(M_PI*M_PI)*FDhalf(mu(V, T)/T)*mu.DV(V, T);
    else FD2V = -std::pow(2.0*mu(V, T), 1.5)/(3.0*M_PI*M_PI)*mu.DV(V, T);

    return FD2V;
}

double FreeEnergy::FDVT(const double V, const double T) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    double FDVleft2  = FDV(V, T - 2*dT);
    double FDVleft1  = FDV(V, T -   dT);
    double FDVright1 = FDV(V, T +   dT);
    double FDVright2 = FDV(V, T + 2*dT);

    return (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dT);
}

double FreeEnergy::FDVT2(const double V, const double T) {

    ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    FermiDirac<ThreeHalf> FD3half;
    FermiDirac<Half>      FDhalf;

    double FDVT;
    double M = mu(V, T);
    double MDT = mu.DT(V, T);

    if (T > 1e-10) 
        FDVT = -std::sqrt(2.0*T)*T/(M_PI*M_PI)*(
                  5.0/3.0*FD3half(M/T)
             + (MDT - M/T)*FDhalf(M/T)
             );
    else FDVT = 0.0;

    return FDVT;
}

double FreeEnergy::FD2T(const double V, const double T) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    double FDTleft2  = FDT(V, T - 2*dT);
    double FDTleft1  = FDT(V, T -   dT);
    double FDTright1 = FDT(V, T +   dT);
    double FDTright2 = FDT(V, T + 2*dT);

    return (-FDTright2 + 8*FDTright1 - 8*FDTleft1 + FDTleft2)/(12.0*dT);
}

double FreeEnergy::FD2T2(const double V, const double T) {

    ChemicalPotential mu;
    mu.setZ(Z);
    mu.setTolerance(tolerance);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    RHSFD2T rhs;
    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSFD2T::dim> yFD2T;

    double xFrom = 1.0;
    double xTo   = 0.0;
    yFD2T.fill(0.0);

    Solver<PD853<RHSFD2T>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    if (T > 1e-10) solver.integrate(rhs, yFD2T, xFrom, xTo);

    return std::pow(Z, -1.0/3.0)*yFD2T[RHSFD2T::result]*rhs.param();
}