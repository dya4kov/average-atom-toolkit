#pragma once
#include <vector>
#include <cstddef>

#include <average-atom-toolkit/configuration.h>

namespace aatk {
namespace atom {

class Atom {
public:
	// new style
	// struct ConstructorArgs {
	// 	double V = 1.0;
	// 	double T = 1.0;
	// 	double Z = 1.0;
	// 	double tolerance = 1e-6;
	// 	int    meshSize = 600;
	// };
	// Atom(ConstructorArgs args);

	// old style
	Atom(
		double V = 1.0, 
		double T = 1.0, 
		double Z = 1.0, 
		double tolerance = 1.e-6,
		int    meshSize = 600
	);
	virtual ~Atom();

	// void reset(ConstructorArgs args);
	virtual void reset(
		double V = 1.0, 
		double T = 1.0, 
		double Z = 1.0, 
		double tolerance = 1.e-6,
		int    meshSize = 600
	);

	double radius();
	double volume();
	double temperature();
	double Znucleus();
	double chemicalPotential();
	double ZfreeElectrons();

	        std::vector<double> U(const std::vector<double>& x);
	        double              U(double x);
	virtual void                U(const double* x, double* y, std::size_t n);

	        std::vector<double> xU(const std::vector<double>& x);
	        double              xU(double x);
	virtual void                xU(const double* x, double* y, std::size_t n);

	        std::vector<double> x2dU(const std::vector<double>& x);
	        double              x2dU(double x);
	virtual void                x2dU(const double* x, double* y, std::size_t n);

	        std::vector<double> electronDensity(const std::vector<double>& x, double eb = 1e+20);
	        double              electronDensity(double x, double eb = 1e+20);
	virtual void                electronDensity(const double* x, double* dens, std::size_t n, double eb = 1e+20);

	        std::vector<double> waveFunction(const std::vector<double>& x, double e, double lambda);
	        double              waveFunction(double x, double e, double lambda);
	virtual void                waveFunction(const double* x, double* y, std::size_t n, double e, double lambda);

	virtual double              innerRP(double e, double lambda);
	virtual double              outerRP(double e, double lambda);
	virtual double              action(double e, double lambda);

	virtual double              energyLevel(int n, int l);

	virtual double              electronStatesDiscrete(int n, int l);
	virtual double              electronStatesDiscrete(int n);
	virtual double              electronStatesDiscrete();

	virtual double              electronStatesContinuous();

protected:

	double V, 
	       T, 
	       Z, 
	       M, 
	       r0, 
	       tolerance;

	int meshSize;

};

}
}