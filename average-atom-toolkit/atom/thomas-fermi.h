#pragma once
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include <average-atom-toolkit/atom/base.h>

namespace aatk {
namespace atom {

class ThomasFermiAtom : public Atom {
public:
	// ThomasFermiAtom(Atom::ConstructorArgs args);

	ThomasFermiAtom(
		double V = 1.0, 
		double T = 1.0, 
		double Z = 1.0, 
		double tolerance = 1.e-6,
		int    meshSize = 600
	);
	~ThomasFermiAtom();

	void reset(
		double V = 1.0, 
		double T = 1.0, 
		double Z = 1.0, 
		double tolerance = 1.e-6,
		int    meshSize = 600
	);

	void U(const double* x, double* y, std::size_t n);
	void xU(const double* x, double* y, std::size_t n);
	void x2dU(const double* x, double* y, std::size_t n);
	void electronDensity(const double* x, double* dens, std::size_t n, double eb = 1e+20);

private:
	std::vector<double> mesh;
	std::vector<double>  pot;
	std::vector<double> dpot;

	gsl_interp_accel  *acc;
	gsl_spline  *phiSpline;
	gsl_spline *dphiSpline;

	void evaluate_chemical_potential();
	void evaluate_potential();

	double mu1(const double V1, const double T1, const double tol);
	double mu1_approx(const double lgV1  /*    T = 0    */);
	double mu1_approx(const double lgV1, const double lgT1);

	// basic buffer
	static std::vector<double>   table;
	static const int    vSize,   tSize;
	static const double lgV0,    lgT0;
	static const double lgVstep, lgTstep;
	static const double bestTolerance;
};

}
}