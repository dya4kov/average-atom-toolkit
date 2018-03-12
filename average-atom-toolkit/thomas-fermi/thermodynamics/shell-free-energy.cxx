#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

#include <average-atom-toolkit/thomas-fermi/thermodynamics/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/ODE/shell-dE.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/shell-chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/shell-free-energy.h>

using numtk::ODE::Array;
using numtk::ODE::Dimension;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using numtk::specfunc::FermiDirac;
using numtk::specfunc::FD::Half;

using aatk::TF::shell::ODE::RHSdE;

using namespace aatk::TF::shell;

FreeEnergy::FreeEnergy() : tolerance(1e-6), Z(1.0), threadsLimit(8) {}

double FreeEnergy::operator() (const double& V, const double& T) {
    return F(V, T);
}

double FreeEnergy::DV(const double& V, const double& T) {
    return FDV(V, T);
}

double FreeEnergy::DT(const double& V, const double& T) {
    return FDT(V, T);
}

double FreeEnergy::D2V(const double& V, const double& T) {
    return FD2V(V, T);
}

double FreeEnergy::DVT(const double& V, const double& T) {
    return FDVT(V, T);
}

double FreeEnergy::D2T(const double& V, const double& T) {
    return FD2T(V, T);
}

// double* FreeEnergy::operator()(
//     const double* V,
//     const double* T, 
//     const std::size_t& vsize, 
//     const std::size_t& tsize
// ) {
//     auto func = std::bind(&FreeEnergy::F, this, _1, _2, _3, _4);
//     return evaluate(func, V, T, vsize, tsize);
// }
// 
// double* FreeEnergy::DV(
//     const double* V,
//     const double* T, 
//     const std::size_t& vsize, 
//     const std::size_t& tsize
// ) {
//     auto func = std::bind(&FreeEnergy::FDV, this, _1, _2, _3, _4);
//     return evaluate(func, V, T, vsize, tsize);
// }
// 
// double* FreeEnergy::DT(
//     const double* V,
//     const double* T, 
//     const std::size_t& vsize, 
//     const std::size_t& tsize
// ) {
//     auto func = std::bind(&FreeEnergy::FDT, this, _1, _2, _3, _4);
//     return evaluate(func, V, T, vsize, tsize);
// }
// 
// double* FreeEnergy::D2V(
//     const double* V, 
//     const double* T, 
//     const std::size_t& vsize, 
//     const std::size_t& tsize
// ) {
//     auto func = std::bind(&FreeEnergy::FD2V, this, _1, _2, _3, _4);
//     return evaluate(func, V, T, vsize, tsize);
// }
// 
// double* FreeEnergy::DVT(
//     const double* V, 
//     const double* T, 
//     const std::size_t& vsize, 
//     const std::size_t& tsize
// ) {
//     auto func = std::bind(&FreeEnergy::FDVT, this, _1, _2, _3, _4);
//     return evaluate(func, V, T, vsize, tsize);
// }
// 
// double* FreeEnergy::D2T(
//     const double* V,
//     const double* T, 
//     const std::size_t& vsize, 
//     const std::size_t& tsize
// ) {
//     auto func = std::bind(&FreeEnergy::FD2T, this, _1, _2, _3, _4);
//     return evaluate(func, V, T, vsize, tsize);
// }
// 
// std::vector<double>& FreeEnergy::operator() (
//     const std::vector<double>& V, 
//     const std::vector<double>& T
// ) {
//     auto Vdata  = V.data();
//     auto Tdata  = T.data();
//     auto Vsize  = V.size();
//     auto Tsize  = T.size();
//     auto func   = std::bind(&FreeEnergy::F, this, _1, _2, _3, _4);
//     auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
//     std::vector<double>* vec_result = 
//          new std::vector<double>(result, result + Vsize*Tsize);
//     return *vec_result;
// }
// 
// std::vector<double>& FreeEnergy::DV(
//     const std::vector<double>& V, 
//     const std::vector<double>& T
// ) {
//     auto Vdata  = V.data();
//     auto Tdata  = T.data();
//     auto Vsize  = V.size();
//     auto Tsize  = T.size();
//     auto func   = std::bind(&FreeEnergy::FDV, this, _1, _2, _3, _4);
//     auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
//     std::vector<double>* vec_result = 
//          new std::vector<double>(result, result + Vsize*Tsize);
//     return *vec_result;
// }
// 
// std::vector<double>& FreeEnergy::DT(
//     const std::vector<double>& V, 
//     const std::vector<double>& T
// ) {
//     auto Vdata  = V.data();
//     auto Tdata  = T.data();
//     auto Vsize  = V.size();
//     auto Tsize  = T.size();
//     auto func   = std::bind(&FreeEnergy::FDT, this, _1, _2, _3, _4);
//     auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
//     std::vector<double>* vec_result = 
//          new std::vector<double>(result, result + Vsize*Tsize);
//     return *vec_result;
// }
// 
// std::vector<double>& FreeEnergy::D2V(
//     const std::vector<double>& V, 
//     const std::vector<double>& T
// ) {
//     auto Vdata  = V.data();
//     auto Tdata  = T.data();
//     auto Vsize  = V.size();
//     auto Tsize  = T.size();
//     auto func   = std::bind(&FreeEnergy::FD2V, this, _1, _2, _3, _4);
//     auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
//     std::vector<double>* vec_result = 
//          new std::vector<double>(result, result + Vsize*Tsize);
//     return *vec_result;
// }
// 
// std::vector<double>& FreeEnergy::DVT(
//     const std::vector<double>& V, 
//     const std::vector<double>& T
// ) {
//     auto Vdata  = V.data();
//     auto Tdata  = T.data();
//     auto Vsize  = V.size();
//     auto Tsize  = T.size();
//     auto func   = std::bind(&FreeEnergy::FDVT, this, _1, _2, _3, _4);
//     auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
//     std::vector<double>* vec_result = 
//          new std::vector<double>(result, result + Vsize*Tsize);
//     return *vec_result;
// }
// 
// std::vector<double>& FreeEnergy::D2T(
//     const std::vector<double>& V, 
//     const std::vector<double>& T
// ) {
//     auto Vdata  = V.data();
//     auto Tdata  = T.data();
//     auto Vsize  = V.size();
//     auto Tsize  = T.size();
//     auto func   = std::bind(&FreeEnergy::FD2T, this, _1, _2, _3, _4);
//     auto result = evaluate(func, Vdata, Tdata, Vsize, Tsize);
//     std::vector<double>* vec_result = 
//          new std::vector<double>(result, result + Vsize*Tsize);
//     return *vec_result;
// }

void FreeEnergy::setZ(const double& _Z) { Z = _Z; }
void FreeEnergy::setThreadsLimit(const std::size_t& Nthreads) {
    threadsLimit = Nthreads;
}
void FreeEnergy::setTolerance(const double& eps) {
    tolerance = eps; 
}

double FreeEnergy::F(const double& V, const double& T) {

    ::aatk::TF::ChemicalPotential M;
    M.setTolerance(tolerance);
    M.setZ(Z);

    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double mu1 = M(V, T)*std::pow(Z, -4.0/3.0);
    
    RHSdE rhs;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSdE::dim> y; y.fill(0.0);

    Solver<PD853<RHSdE>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);
    solver.integrate(rhs, y, 1.0, 0.0);

    double Eint = y[2]*6.0*V1*std::sqrt(2.0) / (M_PI*M_PI)*Z;

    ::aatk::TF::shell::ChemicalPotential dM;
    dM.setTolerance(tolerance);
    dM.setZ(Z);
    dM.setThreadsLimit(threadsLimit);
    return (1.5*Z - Eint)*dM(V, T);
}

double FreeEnergy::FDV(const double& V, const double& T) {

    ::aatk::TF::ChemicalPotential M;
    M.setTolerance(tolerance);
    M.setZ(Z);

    FermiDirac<Half> FDhalf;
    double T1 = T*std::pow(Z, -4.0/3.0);
    double mu = M(V, T);

    double eDens;
    eDens  = T1 > 1e-10 ? 
             std::pow(T, 1.5)*FDhalf(mu/T) :
             std::pow(mu, 1.5)*2.0/3.0;
    eDens *= std::sqrt(2.0) / M_PI / M_PI * Z * Z;

    ::aatk::TF::shell::ChemicalPotential dM;
    dM.setTolerance(tolerance);
    dM.setZ(Z);
    dM.setThreadsLimit(threadsLimit);

    return -eDens*dM(V, T);
}

double FreeEnergy::FDT(const double& V, const double& T) {
    return 0.0;
}

double FreeEnergy::FD2V(const double& V, const double& T) {

    double dV = std::sqrt(std::sqrt(tolerance))*V;

    double FDVleft2  = FDV(V - 2*dV, T);
    double FDVleft1  = FDV(V -   dV, T);
    double FDVright1 = FDV(V +   dV, T);
    double FDVright2 = FDV(V + 2*dV, T);

    return (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dV);
}

double FreeEnergy::FDVT(const double& V, const double& T) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    double FDVleft2  = FDV(V, T - 2*dT);
    double FDVleft1  = FDV(V, T -   dT);
    double FDVright1 = FDV(V, T +   dT);
    double FDVright2 = FDV(V, T + 2*dT);

    return (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dT);
}

double FreeEnergy::FD2T(const double& V, const double& T) {

    double dT = std::sqrt(std::sqrt(tolerance))*T;

    double FDTleft2  = FDT(V, T - 2*dT);
    double FDTleft1  = FDT(V, T -   dT);
    double FDTright1 = FDT(V, T +   dT);
    double FDTright2 = FDT(V, T + 2*dT);

    return (-FDTright2 + 8*FDTright1 - 8*FDTleft1 + FDTleft2)/(12.0*dT);
}