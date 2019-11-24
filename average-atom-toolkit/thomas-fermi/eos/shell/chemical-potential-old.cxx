#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

#include <average-atom-toolkit/thomas-fermi/eos/shell/chemical-potential-old.h>
#include <average-atom-toolkit/thomas-fermi/eos/shell/ODE/dM.h>

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;
using aatk::TF::shell::ODE::RHSdM;

using namespace aatk::TF::shell;

ChemicalPotential::ChemicalPotential() : 
#ifdef ENABLE_MULTITHREADING
    threadsLimit(16), 
#endif
    tolerance(1e-6), Z(1.0)
{}

ChemicalPotential::ChemicalPotential(const ChemicalPotential& dmu) {
    tolerance = dmu.tolerance;
    Z = dmu.Z;
#ifdef ENABLE_MULTITHREADING
    threadsLimit = dmu.threadsLimit;
#endif
}

ChemicalPotential& ChemicalPotential::operator=(const ChemicalPotential& dmu) {
    tolerance = dmu.tolerance;
    Z = dmu.Z;
#ifdef ENABLE_MULTITHREADING
    threadsLimit = dmu.threadsLimit;
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
    tolerance = t;
}

double ChemicalPotential::operator()(const double& V, const double& T) {
    return M(V, T);
}

double* ChemicalPotential::operator()(
    const double* V, 
    const double* T, 
    const std::size_t& vsize, 
    const std::size_t& tsize
) {
    double* result = new double[vsize*tsize];
    for (std::size_t v = 0; v < vsize; ++v) {
        for (std::size_t t = 0; t < tsize; ++t) {
            result[v*tsize + t] = M(V[v], T[t]);
        }
    }
    return result;
}

std::vector<double> ChemicalPotential::operator()(
    const std::vector<double>& V, 
    const std::vector<double>& T 
) {
    std::vector<double> result(V.size()*T.size());
    for (std::size_t v = 0; v < V.size(); ++v) {
        for (std::size_t t = 0; t < T.size(); ++t) {
            result[v*T.size() + t] = M(V[v], T[t]);
        }
    }
    return result;
}

double ChemicalPotential::M(
    const double& V, 
    const double& T
) {
    ::aatk::TF::ElectronStates N;
    N.setVTZ(V, T, Z);
    N.setTolerance(tolerance);
    N.setNmax(25);

#ifdef ENABLE_MULTITHREADING
    N.setThreadsLimit(threadsLimit);
#endif
    
    double eBoundary  = N.eBoundary();
    double continuous = N.continuous(eBoundary);
    double discrete   = N.discrete(eBoundary);
    double dN = discrete - continuous;

    ::aatk::TF::ChemicalPotential M;
    M.setTolerance(tolerance);
    M.setZ(Z);

    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double mu1 = M(V, T)*std::pow(Z, -4.0/3.0);
    
    RHSdM rhs;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    Array<RHSdM::dim> y; y.fill(0.0);

    Solver<PD853<RHSdM>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

	solver.integrate(rhs, y, 1.0, 0.0);
	double dM = y[2]*6.0*V1*std::sqrt(2.0)/(M_PI*M_PI)*std::pow(Z, -1.0/3.0);
	
	return -dN/dM;
}