#pragma once
#include <vector>
#include <cstddef>
#include <average-atom-toolkit/configuration.h>
#ifdef ENABLE_MULTITHREADING
#include <average-atom-toolkit/multithreading/thread-pool.h>
#endif
#include <numeric-toolkit/interpolation/spline.h>

namespace aatk {
namespace semiclassic {

using ::numtk::interpolation::Spline;
#ifdef ENABLE_MULTITHREADING
using ::aatk::multithreading::ThreadPool;
#endif

class Atom {
public:

	Atom(double V = 1.0, 
		 double T = 1.0, 
		 double Z = 1.0, 
		 int    nmax = 10, 
		 bool useContinuous = true,
		 double tolerance = 1.e-6

#ifdef ENABLE_MULTITHREADING
         ,ThreadPool& threads = ::aatk::multithreading::dummy_pool
#endif
	);

	~Atom();

	void                update(double mixing = 0.25);
	void                update(const std::vector<double>& mesh, double mixing = 0.25);
	void                update(const double* mesh, std::size_t size, double mixing = 0.25);
	void                reset(double V = -1.0, double T = -1.0, double Z = -1.0, int nmax = -1);

	double              V();
	double              T();
	double              Z();
	double              M();
	double              N_max();

	std::vector<double> potential(const std::vector<double>& x);
	double              potential(double x);
	void                potential(const double* x, double* y, std::size_t n);

	std::vector<double> waveFunction(double e, double lambda, const std::vector<double>& x);
	double              waveFunction(double e, double lambda, double x);
	void                waveFunction(double e, double lambda, const double* x, double* y, std::size_t n);

	std::vector<double> electronDensity(const std::vector<double>& x);
	double              electronDensity(double x);
	void                electronDensity(const double* x, double* y, std::size_t n);
    double              electronDensityDiscrete(double x);
    double              electronDensityContinuous(double x);

    double              energyLevel(int n, int l);
    double              energyContinuous();
    double              energyFull();

    double              electronStates(int n, int l);
	double              electronStates(int n);
	double              electronStates();
	double              electronStatesContinuous();


	std::array<double, 3> innerRP(double e, double lambda);
	std::array<double, 3> outerRP(double e, double lambda);
	double                action(double e, double lambda);


protected:

	double waveFunctionNorm(
		double energy, double lambda, 
		double rpi,    double rpo, 
		double potrpo, double dpotrpo
	);

	std::vector<double> waveFunctionVec(double e, double lambda, const std::vector<double>& x);

	double electronStates(double chemicalPotential);
    static double electronStatesContinuousFunc (double x, void * classObject);
    double electronStatesContinuous(double chemicalPotential);


    std::vector<double> sorted_mesh(const double* mesh, std::size_t size);
	static double energyDensityContinuousFunc(double x, void * classObject);
	double energyDensityContinuous(double x);
    int evaluateEnergyLevel(int n, int l);
	void evaluateChemicalPotential();

	double volume, temperature, Zcharge, chemPot, r0, nmax, tolerance;

	const std::vector<double>        eLevelStart;
	std::vector<std::vector<double>> eLevel; 
    std::vector<std::vector<bool>>   eLevelReady;
    bool                             chemPotReady;
    bool                             useContinuous;
    int                              nUpdate;

	Spline* densityInterpolation; // sum Nnl |Rnl(x)|^2

#ifdef ENABLE_MULTITHREADING
	ThreadPool& pool;
#endif

};

}
}