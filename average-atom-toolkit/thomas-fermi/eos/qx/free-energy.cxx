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

using namespace aatk::TF::qx;

FreeEnergy::FreeEnergy() : 
#ifdef ENABLE_MULTITHREADING
    threadsLimit(16),
#endif
    tolerance(1e-6),
    Z(1.0), 
    dE0(0.26990017) 
{}

double FreeEnergy::operator()(const double& V, const double& T) {
    double result;
    bool finished;
    F(V, T, result, finished);
    return result;
}

double FreeEnergy::DV(const double& V, const double& T) {
    double result;
    bool finished;
    FDV(V, T, result, finished);
    return result;
}

double FreeEnergy::DT(const double& V, const double& T) {
    double result;
    bool finished;
    FDT(V, T, result, finished);
    return result;
}

double FreeEnergy::D2V(const double& V, const double& T) {
    double result;
    bool finished;
    FD2V(V, T, result, finished);
    return result;
}

double FreeEnergy::DVT(const double& V, const double& T) {
    double result;
    bool finished;
    FDVT(V, T, result, finished);
    return result;
}

double FreeEnergy::D2T(const double& V, const double& T) {
    double result;
    bool finished;
    FD2T(V, T, result, finished);
    return result;
}

double* FreeEnergy::operator()(
    const double* V,
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::F, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::DV(
    const double* V,
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FDV, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::DT(
    const double* V,
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FDT, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::D2V(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FD2V, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::DVT(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FDVT, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::D2T(
    const double* V,
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FD2T, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

std::vector<double> FreeEnergy::operator() (
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func   = std::bind(&FreeEnergy::F, this, _1, _2, _3, _4);
    double* result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    std::vector<double> vresult(result, result + V.size()*T.size());
    delete[] result;
    return vresult;
}

std::vector<double> FreeEnergy::DV(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func   = std::bind(&FreeEnergy::FDV, this, _1, _2, _3, _4);
    double* result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    std::vector<double> vresult(result, result + V.size()*T.size());
    delete[] result;
    return vresult;
}

std::vector<double> FreeEnergy::DT(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func   = std::bind(&FreeEnergy::FDT, this, _1, _2, _3, _4);
    double* result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    std::vector<double> vresult(result, result + V.size()*T.size());
    delete[] result;
    return vresult;
}

std::vector<double> FreeEnergy::D2V(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func   = std::bind(&FreeEnergy::FD2V, this, _1, _2, _3, _4);
    double* result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    std::vector<double> vresult(result, result + V.size()*T.size());
    delete[] result;
    return vresult;
}

std::vector<double> FreeEnergy::DVT(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func   = std::bind(&FreeEnergy::FDVT, this, _1, _2, _3, _4);
    double* result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    std::vector<double> vresult(result, result + V.size()*T.size());
    delete[] result;
    return vresult;
}

std::vector<double> FreeEnergy::D2T(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto func   = std::bind(&FreeEnergy::FD2T, this, _1, _2, _3, _4);
    double* result = evaluate(func, V.data(), T.data(), V.size(), T.size());
    std::vector<double> vresult(result, result + V.size()*T.size());
    delete[] result;
    return vresult;
}

void FreeEnergy::setZ(const double& _Z) { Z = _Z; }

void FreeEnergy::setTolerance(const double& eps) { tolerance = eps; }

#ifdef ENABLE_MULTITHREADING
void FreeEnergy::setThreadsLimit(const std::size_t& N) {
    threadsLimit = std::max(1LU, N);
}

void FreeEnergy::updateThreads(
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

double* FreeEnergy::evaluate(
    std::function<void(const double&, const double&, double&, bool&)> func, 
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    double* result = new double[vsize*tsize];
    bool* finished = new bool[vsize*tsize];
#ifdef ENABLE_MULTITHREADING
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
#else
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            func(V[v], T[t], result[v*tsize + t], finished[v*tsize + t]);
        }
    }
#endif
    delete[] finished;
    return result;
}

void FreeEnergy::F(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {
    
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

    result = dF[RHSdF::result];
    finished = true;
}

void FreeEnergy::FDV(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

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

    if (T > 1e-10) {
        result = -T1*std::sqrt(T1)/(3.0*M_PI*M_PI*M_PI)*(FDhalf(M1/T1)*psi1 + std::sqrt(T1)*Y(M1/T1));
    }
    else {
        result = -std::sqrt(M1)*M1/(9.0*M_PI*M_PI*M_PI)*(2.0*psi1 + 11.0*std::sqrt(M1));
    }

    result *= std::pow(Z, 8.0/3.0); 
    finished = true;
}

void FreeEnergy::FDT(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

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

    result = -dS[RHSdS::result];
    finished = true;
}

void FreeEnergy::FD2V(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

    double dV = std::sqrt(std::sqrt(tolerance))*V;

    bool dummy;

    double FDVleft2 ; FDV(V - 2*dV, T, FDVleft2 , dummy);
    double FDVleft1 ; FDV(V -   dV, T, FDVleft1 , dummy);
    double FDVright1; FDV(V +   dV, T, FDVright1, dummy);
    double FDVright2; FDV(V + 2*dV, T, FDVright2, dummy);

    result = (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dV);
    finished = true;
}

void FreeEnergy::FDVT(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    bool dummy;

    double FDVleft2 ; FDV(V, T - 2*dT, FDVleft2 , dummy);
    double FDVleft1 ; FDV(V, T -   dT, FDVleft1 , dummy);
    double FDVright1; FDV(V, T +   dT, FDVright1, dummy);
    double FDVright2; FDV(V, T + 2*dT, FDVright2, dummy);

    result = (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dT);
    finished = true;
}

void FreeEnergy::FD2T(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    bool dummy;

    double FDTleft2 ; FDT(V, T - 2*dT, FDTleft2 , dummy);
    double FDTleft1 ; FDT(V, T -   dT, FDTleft1 , dummy);
    double FDTright1; FDT(V, T +   dT, FDTright1, dummy);
    double FDTright2; FDT(V, T + 2*dT, FDTright2, dummy);

    result = (-FDTright2 + 8*FDTright1 - 8*FDTleft1 + FDTleft2)/(12.0*dT);
    finished = true;
}