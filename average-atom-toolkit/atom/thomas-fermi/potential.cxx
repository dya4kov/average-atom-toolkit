#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

#include <average-atom-toolkit/util/specfunc/fermi-dirac/complete.h>
#include <average-atom-toolkit/atom/thomas-fermi.h>

namespace aatk {
namespace atom {

using ::aatk::util::specfunc::FermiDirac;
using ::aatk::util::specfunc::FD::Half;

static const int TFpotentialODEdim = 2;

struct TFpotentialODEParams {
	double mu;
	double a;
	double T;
	FermiDirac<Half> FDhalf;
};

int TFpotentialRHS_T(double x, const double y[], double dydx[], void *params) {
	auto& p      = *((TFpotentialODEParams *) params);
	auto& mu     = p.mu;
	auto& a      = p.a;
	auto& T      = p.T;
	auto& FDhalf = p.FDhalf;

	dydx[0] = 2.0*x*y[1];
	double phi = y[0] + x*x*mu;
	if (x > 0) {
		dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf(phi/(T*x*x));
	}
	else {
		dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
	}

	return GSL_SUCCESS;
}

int TFpotentialRHS_0(double x, const double y[], double dydx[], void *params) {
	auto& p  = *((TFpotentialODEParams *) params);
	auto& mu = p.mu;
	auto& a  = p.a;

	double phi = y[0] + x*x*mu;
	dydx[0] = 2.0*x*y[1];
	dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
	return GSL_SUCCESS;
}

void integrate(
	int   (*rhs)(double, const double[], double[], void*), 
	double *phi,
	double  xFrom, 
	double  xTo,
	double  tolerance,
	void   *params
) {

	gsl_odeiv2_step*    stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, TFpotentialODEdim);
	gsl_odeiv2_control* control = gsl_odeiv2_control_y_new (0.0, 0.1*tolerance);
	gsl_odeiv2_evolve*   evolve = gsl_odeiv2_evolve_alloc (TFpotentialODEdim);

	gsl_odeiv2_system sys = {rhs, NULL, TFpotentialODEdim, params};

	double h = -std::min(1.e-6, tolerance);
	while(xFrom > xTo) {
		int status = gsl_odeiv2_evolve_apply (
			evolve, control, stepper,
		    &sys,
		    &xFrom, xTo,
		    &h, phi
		);
		if (status != GSL_SUCCESS) {
		    std::cout << "ODE for potential integration problems" << std::endl;
		}
	}

	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (stepper);
}

void ThomasFermiAtom::evaluate_potential() {
	const double* u = mesh.data();
	double*       y = pot.data();
	double*      dy = dpot.data();

	// u - sorted array from 0.0 to 1.0
	double mu1 = M*std::pow(Z, -4.0/3.0);
	double  V1 = V*Z;
	double  T1 = T*std::pow(Z, -4.0/3.0);

	TFpotentialODEParams params;
	params.mu = mu1;
	params.T  = T1;
	params.a  = std::pow(2.0,   7.0/6.0)
	          * std::pow(3.0,   2.0/3.0)
	          * std::pow(M_PI, -5.0/3.0)
	          * std::pow(V1,    2.0/3.0);

	double phi[2] = { 0 }; 

	auto rhs = T1 <= 1e-10 ? &TFpotentialRHS_0 : &TFpotentialRHS_T;

	gsl_odeiv2_step*    stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, TFpotentialODEdim);
	gsl_odeiv2_control* control = gsl_odeiv2_control_y_new (0.0, 0.1*tolerance);
	gsl_odeiv2_evolve*   evolve = gsl_odeiv2_evolve_alloc (TFpotentialODEdim);

	gsl_odeiv2_system sys = {rhs, NULL, TFpotentialODEdim, &params};

	double from = 1.0;
	double h = -std::min(1.e-6, tolerance);

	for (int i = meshSize - 1; i >= 0; --i) {
		double to = u[i];
		while(from > to) {
			int status = gsl_odeiv2_evolve_apply (
				evolve, control, stepper,
				&sys,
				&from, to,
				&h, phi
			);
			if (status != GSL_SUCCESS) {
				std::cerr << "ODE for potential integration problems" << std::endl;
			}
		}
		y[i] = phi[0];
		dy[i] = phi[1];
		from = to;
	}

	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (stepper);

	// potential and its derivative interpolation
	if (phiSpline != nullptr) gsl_spline_free(phiSpline);
	if (dphiSpline != nullptr) gsl_spline_free(dphiSpline);
    phiSpline = gsl_spline_alloc(gsl_interp_cspline, meshSize);
    dphiSpline = gsl_spline_alloc(gsl_interp_cspline, meshSize);
    gsl_spline_init(phiSpline, u, y, meshSize);
    gsl_spline_init(dphiSpline, u, dy, meshSize);

	return;
}

// init static variables
const int    ThomasFermiAtom::vSize          =   181;
const int    ThomasFermiAtom::tSize          =   201;
const double ThomasFermiAtom::lgV0           = -10.0;
const double ThomasFermiAtom::lgT0           = -10.0;
const double ThomasFermiAtom::lgVstep        =   0.1;
const double ThomasFermiAtom::lgTstep        =   0.1;
const double ThomasFermiAtom::bestTolerance  = 1e-12;

double ThomasFermiAtom::mu1(const double V1, const double T1, const double tol) {

	TFpotentialODEParams params;
	params.T = T1;
	params.a = std::pow(2.0, 7.0/6.0)
	         * std::pow(3.0, 2.0/3.0)
	         * std::pow(M_PI, -5.0/3.0)
	         * std::pow(V1, 2.0/3.0);

	double phi_0 = std::pow(4.0*M_PI/3.0/V1, 1.0/3.0);

	auto rhs = T1 <= 1e-10 ? &TFpotentialRHS_0 : &TFpotentialRHS_T;

	double muStart;
	if (T1 <= 1e-10) muStart = mu1_approx(std::log10(V1));
	else muStart = mu1_approx(std::log10(V1), std::log10(T1));

	double delta = 0.1;

	double muPrev = muStart;
	params.mu = muPrev;

	double phi[2] = {0.0};
	integrate(rhs, phi, 1.0, 0.0, tolerance, (void*) &params);
	double phiPrev = phi[0];

	double muCurr = muStart - delta*std::abs(muStart) - tolerance;
	params.mu = muCurr;

	phi[0] = phi[1] = 0.0;
	integrate(rhs, phi, 1.0, 0.0, tolerance, (void*) &params);
	double phiCurr = phi[0];

	double muNext = 0.0;
	double phiNext = 0.0;

	double error = std::abs(muPrev - muCurr)/std::abs(muPrev + muCurr + tolerance);

	while (error > tolerance) {
		muNext = muCurr - (phiCurr - phi_0)*(muCurr - muPrev)/(phiCurr - phiPrev);
		params.mu = muNext;

		phi[0] = phi[1] = 0.0;
		integrate(rhs, phi, 1.0, 0.0, tolerance, (void*) &params);
		phiNext = phi[0];

		muPrev  = muCurr;
		phiPrev = phiCurr;
		muCurr  = muNext;
		phiCurr = phiNext;

		error = std::abs(muPrev - muCurr)/std::abs(muPrev + muCurr + tolerance);
	}

	return muNext;
}

double ThomasFermiAtom::mu1_approx(const double lgV) {
	const double B0 =  0.648742997083556;
	const double B1 = -0.704984628856768;
	const double B2 = -0.0224226496439102;
	const double B3 = -0.00419385235723519;
	const double B4 = -3.75915351702641E-5;
	const double B5 =  3.94764845762704E-5;
	const double B6 =  5.4828018180471E-7;
	const double B7 = -1.49964096611993E-7;
	double lgV1 = lgV;
	double lgV2 = lgV1*lgV1;
	double lgV3 = lgV1*lgV2;
	double lgV4 = lgV2*lgV2;
	double lgV5 = lgV3*lgV2;
	double lgV6 = lgV3*lgV3;
	double lgV7 = lgV4*lgV3;
	double lgMu = B0 + B1*lgV1 + B2*lgV2 + B3*lgV3 + B4*lgV4 + B5*lgV5 + B6*lgV6 + B7*lgV7;
	return std::pow(10.0, lgMu);
}

double ThomasFermiAtom::mu1_approx(const double lgV, const double lgT) {
	int v, t;
	v = (int) std::floor((lgV - lgV0)/lgVstep);
	t = (int) std::floor((lgT - lgT0)/lgTstep);
	double result = 0.0;
	if ((v >= 0 && v < vSize) || (t < tSize)) {
		double f00, f01, f10, f11;
		f00 = table[v     +       t*vSize];
		f01 = table[v     + (t + 1)*vSize];
		f10 = table[v + 1 +       t*vSize];
		f11 = table[v + 1 + (t + 1)*vSize];
		int sign = 0;
		sign = (sign == 0 && f00 > 0.0 && f01 > 0.0 && f10 > 0.0 && f11 > 0.0) ?  1 : 0;
		sign = (sign == 0 && f00 < 0.0 && f01 < 0.0 && f10 < 0.0 && f11 < 0.0) ? -1 : 0;
		if (sign != 0) {
			f00 = std::log10(std::abs(f00));
			f01 = std::log10(std::abs(f01));
			f10 = std::log10(std::abs(f10));
			f11 = std::log10(std::abs(f11));
		}
		double V0 = (lgV - lgV0)/lgVstep - 1.0*v;
		double T0 = (lgT - lgT0)/lgTstep - 1.0*t;
		double a[2][2];
		a[0][0] = f00;
		a[1][0] = f10 - f00;
		a[0][1] = f01 - f00;
		a[1][1] = f11 + f00 - (f10 + f01);
		for (int i = 0; i < 2; ++i) {
		    for (int j = 0; j < 2; ++j) {
		        result += a[i][j]*std::pow(V0, i)*std::pow(T0, j);
		    }
		}
		if (sign != 0) result = sign*std::pow(10.0, result);
	}
	return result;
}

}
}