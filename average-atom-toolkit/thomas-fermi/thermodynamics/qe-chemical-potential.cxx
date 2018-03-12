#include <cmath>
#include <thread>

#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

#include <average-atom-toolkit/thomas-fermi/quantum-exchange/potential.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/qe-chemical-potential.h>

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::MHalf;

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using std::placeholders::_4;

using namespace aatk::TF::QE;

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

void ChemicalPotential::M(
    const double& V, 
    const double& T, 
    double& result, 
    bool& finished
) {

    ::aatk::TF::QE::Potential    psi;
    psi.setVTZ(V, T, Z);
    psi.setTolerance(tolerance);

    ::aatk::TF::ChemicalPotential mu;
    mu.setZ(Z);    
    mu.setTolerance(tolerance);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    FermiDirac<MHalf> FDmhalf;

    double dM = std::sqrt(2.0)/ (6.0 * M_PI);

    if (T > 1e-10) {
        dM *= 0.5*std::sqrt(T)*FDmhalf(mu1/T1) + psi(1.0);
    }
    else {
        dM *= 1.0 + psi(1.0);
    }

    result = dM*std::pow(Z, 2.0/3.0);
    finished = true;
}