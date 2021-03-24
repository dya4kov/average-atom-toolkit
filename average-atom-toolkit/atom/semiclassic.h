#pragma once
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include <average-atom-toolkit/atom/base.h>

namespace aatk {
namespace atom {

class SemiclassicAtom : public Atom {
public:
	SemiclassicAtom(
		double V = 1.0, 
		double T = 1.0, 
		double Z = 1.0, 
		double tolerance = 1.e-6,
		int    meshSize = 1600,
		int    nmax = 20,
		bool   useContinuous = true
	);
	~SemiclassicAtom();

	void reset(
		double V = 1.0, 
		double T = 1.0, 
		double Z = 1.0, 
		double tolerance = 1.e-6,
		int    meshSize = 1600,
		int    nmax = 20,
		bool   useContinuous = true
	);

	void   update(double mixing = 0.75);

	double U(double x);
	double electronDensity(double x);

	void   U(const double* x, double* y, std::size_t n);
	void   xU(const double* x, double* y, std::size_t n);
	void   x2dU(const double* x, double* y, std::size_t n);
	void   electronDensity(const double* x, double* dens, std::size_t n, double eb = 1e+20);
	void   waveFunction(const double* x, double* y, std::size_t n, double e, double lambda);

	double innerRP(double e, double lambda);
	double outerRP(double e, double lambda);
	double action(double e, double lambda);

	double energyLevel(int n, int l);

	std::vector<double> electronDensityContinuous(const std::vector<double>& x);
    double              electronDensityContinuous(double x);
    void                electronDensityContinuous(const double* x, double* y, std::size_t n);


	double electronStatesDiscrete(int n, int l);
	double electronStatesDiscrete(int n);
	double electronStatesDiscrete();
	double electronStatesContinuous();

	double electronStatesContinuous(double CP);
	double electronStatesDiscrete(double chemicalPotential);

	double boundaryEnergyValue();

	int discreteLevelsNumber();

	double energyFull();
	double energyContinuous();
	double internalEnergy();

	double entropy();

private:

	bool useContinuous;
	int nmax, nUpdate;
	double boundaryEnergy;
	std::vector<double> mesh;
	std::vector<double>  pot;
	std::vector<double> dpot;
	std::vector<double> dens;
	
	const std::vector<double>        eLevelStart;
	std::vector<std::vector<double>> eLevel; 
    std::vector<std::vector<bool>>    eLevelReady; // was int !!
    bool                             chemPotReady;

	gsl_interp_accel *phiAcc;
	gsl_interp_accel *dphiAcc;
	gsl_interp_accel *densAcc;

	gsl_spline  *phiSpline;
	gsl_spline *dphiSpline;
	gsl_spline *densSpline;

	void evaluate_potential();
	void evaluate_chemical_potential();
	void evaluate_energy_level(int n, int l);
    void evaluate_boundary_energy();
	
	//?double electronStatesContinuousFunc (double x, void * atomClass);


	void evaluate_wave_function(
		const double* u, 
		double*       wf, 
		std::size_t   n,
		double        energy, 
		double        lambda
	);

	void normalize_wave_function(
		const double* u, 
		double*       wf, 
		std::size_t   n
	);

};

}
}