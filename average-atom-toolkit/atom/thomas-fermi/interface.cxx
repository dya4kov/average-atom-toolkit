#include <cmath>
#include <iostream>
#include <average-atom-toolkit/util/specfunc/fermi-dirac/complete.h>
#include <average-atom-toolkit/util/specfunc/fermi-dirac/incomplete.h>

#include <average-atom-toolkit/atom/thomas-fermi.h>

namespace aatk {
namespace atom {

using aatk::util::specfunc::FermiDirac;
using aatk::util::specfunc::FermiDiracInc;
using aatk::util::specfunc::FD::Half;
using aatk::util::specfunc::FDI::HalfInc;

// ThomasFermiAtom::ThomasFermiAtom(Atom::ConstructorArgs args) :
// 	Atom(args), 
// 	acc(nullptr),
//     phiSpline(nullptr),
//     dphiSpline(nullptr)
// {
// 	acc = gsl_interp_accel_alloc();
// 	reset(args);
// }

// old style
ThomasFermiAtom::ThomasFermiAtom(
	double _V,
	double _T,
	double _Z,
	double _tolerance,
	int    _meshSize
) : 
	Atom(_V, _T, _Z, _tolerance, _meshSize),
	acc(nullptr),
	phiSpline(nullptr),
	dphiSpline(nullptr)
{
	acc = gsl_interp_accel_alloc();
	reset(_V, _T, _Z, _tolerance, _meshSize);
}

void ThomasFermiAtom::reset(
	double _V,
	double _T,
	double _Z,
	double _tolerance,
	int    _meshSize
) {
	V = _V;
	T = _T;
	Z = _Z;
	tolerance = _tolerance;
	meshSize = _meshSize;
	r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
	evaluate_chemical_potential();
	// setup mesh
	mesh.resize(meshSize);
	double umin = 0.0;
	double umax = 1.0;
	for (int i = 0; i < meshSize; ++i) {
		mesh[i] = umin + i*(umax - umin)/(meshSize - 1);
	}
	pot.resize(meshSize, 0.0);
	dpot.resize(meshSize, 0.0);
	evaluate_potential();
}

ThomasFermiAtom::~ThomasFermiAtom() {
	gsl_spline_free(phiSpline);
	gsl_spline_free(dphiSpline);
	gsl_interp_accel_free(acc);
}

void ThomasFermiAtom::evaluate_chemical_potential() { 
	double V1 = V*Z;
	double T1 = T*std::pow(Z, -4.0/3.0);
	double M1 = mu1(V1, T1, tolerance);
	M = M1*std::pow(Z, 4.0/3.0);
}

void ThomasFermiAtom::U(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i]  = gsl_spline_eval(phiSpline, std::sqrt(x[i]), acc);
		y[i] *= std::pow(Z, 4.0/3.0)/x[i];
	}
	return;
}

void ThomasFermiAtom::xU(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i]  = gsl_spline_eval(phiSpline, std::sqrt(x[i]), acc);
		y[i] *= std::pow(Z, 4.0/3.0);
	}
	return;
}

void ThomasFermiAtom::x2dU(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		double u = std::sqrt(x[i]);
		double y0 = gsl_spline_eval(phiSpline, u, acc);
		double y1 = gsl_spline_eval(dphiSpline, u, acc);
		y[i]  = 0.5*y1*u + y0;
		y[i] *= std::pow(Z, 4.0/3.0);
	}
	return;
}

void ThomasFermiAtom::electronDensity(const double* x, double* dens, std::size_t n, double eb) {

	double Zm43 = std::pow(Z, -4.0/3.0);

	double T1  = T*Zm43;
	double M1  = M*Zm43;
	double eb1 = eb*Zm43;

	FermiDirac<Half> FDhalf;
	FermiDiracInc<HalfInc> FDhalfI;

	for (std::size_t i = 0; i < n; ++i) {
		double ypot  = gsl_spline_eval(phiSpline, std::sqrt(x[i]), acc);
		double phiMu = ypot + x[i]*M1;
		double phiEb = ypot + x[i]*eb1;

		double result;

		if (phiEb <= 0.0) result = 0.0;
		else {
			double argX = phiMu/(T1*x[i]);
			double argY = phiEb/(T1*x[i]);
			if (argY - argX > 25.0)
				result = T*std::sqrt(2.0*T)*FDhalf(argX)/(M_PI*M_PI);
			else
				result = T*std::sqrt(2.0*T)*FDhalfI(argX, argY)/(M_PI*M_PI);
		}
		dens[i] = result;
	}
	return;
}

}
}