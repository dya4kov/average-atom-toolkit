#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/thermodynamics/ODE/shell-dM.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/shell-chemical-potential.h>

using namespace aatk::TF::shell;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;
using aatk::TF::shell::ODE::RHSdM;

ChemicalPotential::ChemicalPotential() : 
    tolerance(1e-6), Z(1.0), threadsLimit(4) {}

ChemicalPotential::ChemicalPotential(const ChemicalPotential& dmu) {
    tolerance = dmu.tolerance;
    Z = dmu.Z;
    threadsLimit = dmu.threadsLimit;
    mu = dmu.mu;
}

ChemicalPotential& ChemicalPotential::operator=(const ChemicalPotential& dmu) {
    tolerance = dmu.tolerance;
    threadsLimit = dmu.threadsLimit;
    Z = dmu.Z;
    mu = dmu.mu;
    return *this;
}

void ChemicalPotential::setTolerance(const double& t) {
    tolerance = t;
    mu.setTolerance(t);
    N .setTolerance(t);
}

void ChemicalPotential::setThreadsLimit(const std::size_t& Nthreads) {
    threadsLimit = Nthreads; 
}

double ChemicalPotential::operator()(const double& V, const double& T) {
    double result;
    bool finished;
    M(V, T, result, finished);
    return result;
}

double ChemicalPotential::M(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

    N.set(V); N.set(T);
 
    double eBoundary  = N.eBoundary();
    double continuous = N.continuous(eBoundary);
    double discrete   = N.discrete(eBoundary);
    double dN = discrete - continuous;

    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    
    RHSdM rhs;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSdM::dim> y; y.fill(0.0);

    Solver<PD853<RHSdM>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

	solver.integrate(rhs, y, 1.0, 0.0);
	double dM = y[2]*6.0*V1*std::sqrt(2.0)/(M_PI*M_PI)*pow(Z, -1.0/3.0);
	
	result = -dN/dM;
	finished = true;
}