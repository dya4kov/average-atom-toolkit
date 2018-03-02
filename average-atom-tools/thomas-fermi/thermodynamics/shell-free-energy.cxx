#include <cmath>
#include <thread>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

#include <average-atom-tools/thomas-fermi/thermodynamics/chemical-potential.h>
#include <average-atom-tools/thomas-fermi/thermodynamics/shell-chemical-potential.h>
#include <average-atom-tools/thomas-fermi/thermodynamics/ODE/shell-dE.h>

using numtools::ODE::Array;
using numtools::ODE::Dimension;
using numtools::ODE::Solver;
using numtools::ODE::stepper::PD853;

using numtools::specfunc::FermiDirac;
using numtools::specfunc::FD::Half;

using AATools::TF::shell::ODE::RHSdE;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;

using namespace AATools::TF::shell;

FreeEnergy::FreeEnergy() : tolerance(1e-6), Z(1.0), threadsLimit(4) {}

double FreeEnergy::operator() (const double& V, const double& T) {

    ::AATools::TF::ChemicalPotential         M;
    ::AATools::TF::shell::ChemicalPotential dM;

     M.setTolerance(tolerance);
    dM.setTolerance(tolerance);

     M.setZ(Z);
    dM.setZ(Z);

     M.setThreadsLimit(threadsLimit);
    dM.setThreadsLimit(threadsLimit);

    double result;
    bool finished;
    F(V, T, result, finished);
    return result;
}

double FreeEnergy::DV(const double& V, const double& T) {

    ::AATools::TF::ChemicalPotential         M;
    ::AATools::TF::shell::ChemicalPotential dM;

     M.setTolerance(tolerance);
    dM.setTolerance(tolerance);

     M.setZ(Z);
    dM.setZ(Z);

     M.setThreadsLimit(threadsLimit);
    dM.setThreadsLimit(threadsLimit);

    double result;
    bool finished;
    FDV(V, T, result, finished);
    return result;
}

double FreeEnergy::DT(const double& V, const double& T) {

    ::AATools::TF::ChemicalPotential         M;
    ::AATools::TF::shell::ChemicalPotential dM;

     M.setTolerance(tolerance);
    dM.setTolerance(tolerance);

     M.setZ(Z);
    dM.setZ(Z);

     M.setThreadsLimit(threadsLimit);
    dM.setThreadsLimit(threadsLimit);

    double result;
    bool finished;
    FDT(V, T, result, finished);
    return result;
}

double FreeEnergy::D2V(const double& V, const double& T) {

    ::AATools::TF::ChemicalPotential         M;
    ::AATools::TF::shell::ChemicalPotential dM;

     M.setTolerance(tolerance);
    dM.setTolerance(tolerance);

     M.setZ(Z);
    dM.setZ(Z);

     M.setThreadsLimit(threadsLimit);
    dM.setThreadsLimit(threadsLimit);

    double result;
    bool finished;
    FD2V(V, T, result, finished);
    return result;
}

double FreeEnergy::DVT(const double& V, const double& T) {

    ::AATools::TF::ChemicalPotential         M;
    ::AATools::TF::shell::ChemicalPotential dM;

     M.setTolerance(tolerance);
    dM.setTolerance(tolerance);

     M.setZ(Z);
    dM.setZ(Z);

     M.setThreadsLimit(threadsLimit);
    dM.setThreadsLimit(threadsLimit);
    
    double result;
    bool finished;
    FDVT(V, T, result, finished);
    return result;
}

double FreeEnergy::D2T(const double& V, const double& T) {

    ::AATools::TF::ChemicalPotential         M;
    ::AATools::TF::shell::ChemicalPotential dM;

     M.setTolerance(tolerance);
    dM.setTolerance(tolerance);

     M.setZ(Z);
    dM.setZ(Z);

     M.setThreadsLimit(threadsLimit);
    dM.setThreadsLimit(threadsLimit);

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

std::vector<double>& FreeEnergy::operator() (
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto Vdata  = V.data();
    auto Tdata  = T.data();
    auto Vsize  = V.size();
    auto Tsize  = T.size();
    auto func   = std::bind(&FreeEnergy::F, this, _1, _2, _3, _4);
    auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
    std::vector<double>* vec_result = 
         new std::vector<double>(result, result + Vsize*Tsize);
    return *vec_result;
}

std::vector<double>& FreeEnergy::DV(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto Vdata  = V.data();
    auto Tdata  = T.data();
    auto Vsize  = V.size();
    auto Tsize  = T.size();
    auto func   = std::bind(&FreeEnergy::FDV, this, _1, _2, _3, _4);
    auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
    std::vector<double>* vec_result = 
         new std::vector<double>(result, result + Vsize*Tsize);
    return *vec_result;
}

std::vector<double>& FreeEnergy::DT(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto Vdata  = V.data();
    auto Tdata  = T.data();
    auto Vsize  = V.size();
    auto Tsize  = T.size();
    auto func   = std::bind(&FreeEnergy::FDT, this, _1, _2, _3, _4);
    auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
    std::vector<double>* vec_result = 
         new std::vector<double>(result, result + Vsize*Tsize);
    return *vec_result;
}

std::vector<double>& FreeEnergy::D2V(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto Vdata  = V.data();
    auto Tdata  = T.data();
    auto Vsize  = V.size();
    auto Tsize  = T.size();
    auto func   = std::bind(&FreeEnergy::FD2V, this, _1, _2, _3, _4);
    auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
    std::vector<double>* vec_result = 
         new std::vector<double>(result, result + Vsize*Tsize);
    return *vec_result;
}

std::vector<double>& FreeEnergy::DVT(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto Vdata  = V.data();
    auto Tdata  = T.data();
    auto Vsize  = V.size();
    auto Tsize  = T.size();
    auto func   = std::bind(&FreeEnergy::FDVT, this, _1, _2, _3, _4);
    auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
    std::vector<double>* vec_result = 
         new std::vector<double>(result, result + Vsize*Tsize);
    return *vec_result;
}

std::vector<double>& FreeEnergy::D2T(
    const std::vector<double>& V, 
    const std::vector<double>& T
) {
    auto Vdata  = V.data();
    auto Tdata  = T.data();
    auto Vsize  = V.size();
    auto Tsize  = T.size();
    auto func   = std::bind(&FreeEnergy::FD2T, this, _1, _2, _3, _4);
    auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
    std::vector<double>* vec_result = 
         new std::vector<double>(result, result + Vsize*Tsize);
    return *vec_result;
}

void FreeEnergy::setZ(const double& _Z) { Z = _Z; }
void FreeEnergy::setThreadsLimit(const std::size_t& N) {
    threadsLimit = std::max(1LU, N);
}
void FreeEnergy::setTolerance(const double& eps) { tolerance = eps; }

double* FreeEnergy::evaluate(
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
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            finished[v*tsize + t] = false;
            std::thread run(func, std::cref(V[v]), std::cref(T[t]), 
                                  std::ref(result[v*tsize + t]), 
                                  std::ref(finished[v*tsize + t]));
            run.detach(); ++threads;
            while (threads == threadsLimit) {
                while (finished[current]) {
                    --threads;
                    ++current;
                }
            }
        }
    }
    bool all_finished = false;
    while (!all_finished) {
        while (finished[current]) {
            ++current; if (current == vsize*tsize) break;
        }
        all_finished = (current == vsize*tsize);
    }
    delete[] finished;
    return result;
}

void FreeEnergy::F(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    
    RHSdE rhs;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSdE::dim> y; y.fill(0.0);

    Solver<PD853<RHSdE>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    solver.integrate(rhs, y, 1.0, 0.0);

    double Eint = y[2]*6.0*V1*std::sqrt(2.0) / (M_PI*M_PI)*Z;

    ChemicalPotential dM;
    result = (1.5*Z - Eint)*dM(V, T);
    finished = true;
}

void FreeEnergy::FDV(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished) 
{

    ::AATools::TF::ChemicalPotential         M;
    ::AATools::TF::shell::ChemicalPotential dM;

    FermiDirac<Half> FDhalf;
    double T1 = T*std::pow(Z, -4.0/3.0);

    double eDens;
    eDens  = T1 > 1e-10 ? 
             std::pow(T, 1.5)*FDhalf(M(V, T)/T) :
             std::pow(M(V, T), 1.5)*2.0/3.0;
    eDens *= std::sqrt(2.0) / M_PI / M_PI * Z * Z;

    result = -eDens*dM(V, T);
    finished = true;
}

void FreeEnergy::FDT(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {
    result = 0.0;
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