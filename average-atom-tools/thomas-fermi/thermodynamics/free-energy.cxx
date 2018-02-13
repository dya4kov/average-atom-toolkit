#include <cmath>
#include <thread>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

#include <average-atom-tools/thomas-fermi/atom/potential.h>
#include <average-atom-tools/thomas-fermi/thermodynamics/free-energy.h>
#include <average-atom-tools/thomas-fermi/thermodynamics/ODE/F.h>
#include <average-atom-tools/thomas-fermi/thermodynamics/ODE/FDT.h>

using numtools::ODE::Array;
using numtools::ODE::Dimension;
using numtools::ODE::Solver;
using numtools::ODE::stepper::PD853;

using ::numtools::specfunc::FermiDirac;
using ::numtools::specfunc::FD::ThreeHalf;

using AATools::TF::ODE::RHSF1;
using AATools::TF::ODE::RHSF2;
using AATools::TF::ODE::RHSFDT;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;

using namespace AATools::TF;

FreeEnergy::FreeEnergy() : tolerance(1e-6), Z(1.0), E0(0.76874512422), threadsLimit(4) {}

double FreeEnergy::operator() (const double& V, const double& T) {
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
	const size_t& vsize, 
	const size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::F, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::DV(
	const double* V,
	const double* T, 
	const size_t& vsize, 
	const size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FDV, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::DT(
	const double* V,
	const double* T, 
	const size_t& vsize, 
	const size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FDT, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::D2V(
	const double* V, 
	const double* T, 
	const size_t& vsize, 
	const size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FD2V, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::DVT(
	const double* V, 
	const double* T, 
	const size_t& vsize, 
	const size_t& tsize
) {
    auto func = std::bind(&FreeEnergy::FDVT, this, _1, _2, _3, _4);
    return evaluate(func, V, T, vsize, tsize);
}

double* FreeEnergy::D2T(
	const double* V,
	const double* T, 
	const size_t& vsize, 
	const size_t& tsize
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
void FreeEnergy::setThreadsLimit(const size_t& N) {
    threadsLimit = std::max(1LU, N);
}
void FreeEnergy::setTolerance(const double& eps) { tolerance = eps; }

double* FreeEnergy::evaluate(
    std::function<void(const double&, const double&, double&, bool&)> func, 
    const double* V, 
    const double* T, 
    const size_t& vsize, 
    const size_t& tsize
) {
    double* result = new double[vsize*tsize];
    bool* finished = new bool[vsize*tsize];
    size_t threads = 0;
    size_t current = 0;
    for (size_t v = 0; v < vsize; ++v) {
        for (size_t t = 0; t < tsize; ++t) {
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

void FreeEnergy::F(const double& V, const double& T, double& result, bool& finished) {

    Potential phi;
    phi.setV(V); 
    phi.setT(T);
    phi.setZ(Z);
    phi.setTolerance(tolerance);

    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double muZ = phi.mu();
    double mu1 = muZ*std::pow(Z, -4.0/3.0);

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

    double F = E0 + F1[RHSF1::result]*rhs1.param()
                  + F2[RHSF2::result]*rhs2.param();

    F *= std::pow(Z, 7.0/3.0);
    F += 0.5*muZ*Z;

    result = F; finished = true;
}

void FreeEnergy::FDV(const double& V, const double& T, double& result, bool& finished) {

    Potential phi;
    phi.setV(V); 
    phi.setT(T);
    phi.setZ(Z);
    phi.setTolerance(tolerance);

    double mu = phi.mu();

    FermiDirac<ThreeHalf> FD3half;

    double FDV;
    if (T > 1e-10) FDV = -std::pow(2.0*T, 2.5)/(6.0*M_PI*M_PI)*FD3half(mu/T);
    else FDV = -std::pow(2.0*mu, 2.5)/(15.0*M_PI*M_PI);

    result = FDV; finished = true;
}

void FreeEnergy::FDT(const double& V, const double& T, double& result, bool& finished) {

	Potential phi;
    phi.setV(V); 
    phi.setT(T);
    phi.setZ(Z);
    phi.setTolerance(tolerance);

    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double muZ = phi.mu();
    double mu1 = muZ*std::pow(Z, -4.0/3.0);

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

	result = FDT[RHSFDT::result]*rhs.param();
	finished = true;
}

void FreeEnergy::FD2V(const double& V, const double& T, double& result, bool& finished) {

	double dV = std::sqrt(std::sqrt(tolerance))*V;

	bool dummy;

	double FDVleft2 ; FDV(V - 2*dV, T, FDVleft2 , dummy);
	double FDVleft1 ; FDV(V -   dV, T, FDVleft1 , dummy);
	double FDVright1; FDV(V +   dV, T, FDVright1, dummy);
	double FDVright2; FDV(V + 2*dV, T, FDVright2, dummy);

	result = (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dV);
	finished = true;
}

void FreeEnergy::FDVT(const double& V, const double& T, double& result, bool& finished) {

	double dT = std::sqrt(std::sqrt(tolerance))*T;

	bool dummy;

	double FDVleft2 ; FDV(V, T - 2*dT, FDVleft2 , dummy);
	double FDVleft1 ; FDV(V, T -   dT, FDVleft1 , dummy);
	double FDVright1; FDV(V, T +   dT, FDVright1, dummy);
	double FDVright2; FDV(V, T + 2*dT, FDVright2, dummy);

	result = (-FDVright2 + 8*FDVright1 - 8*FDVleft1 + FDVleft2)/(12.0*dT);
	finished = true;
}

void FreeEnergy::FD2T(const double& V, const double& T, double& result, bool& finished) {

	double dT = std::sqrt(std::sqrt(tolerance))*T;

	bool dummy;

	double FDTleft2 ; FDT(V, T - 2*dT, FDTleft2 , dummy);
	double FDTleft1 ; FDT(V, T -   dT, FDTleft1 , dummy);
	double FDTright1; FDT(V, T +   dT, FDTright1, dummy);
	double FDTright2; FDT(V, T + 2*dT, FDTright2, dummy);

	result = (-FDTright2 + 8*FDTright1 - 8*FDTleft1 + FDTleft2)/(12.0*dT);
	finished = true;
}