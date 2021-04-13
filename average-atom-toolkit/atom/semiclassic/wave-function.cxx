#include <cmath>
#include <algorithm>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

#include <average-atom-toolkit/util/specfunc/bessel/Jnu.h>
#include <average-atom-toolkit/util/specfunc/bessel/Ynu.h>
#include <average-atom-toolkit/util/specfunc/bessel/Knu.h>

#include <average-atom-toolkit/atom/semiclassic.h>

namespace aatk {
namespace atom {

void SemiclassicAtom::waveFunction(
	const double *x,
	double       *result,
	std::size_t   n,
	double        energy, 
	double        lambda
) {
	std::vector<double> u(n);
	for (int i = 0; i < u.size(); ++i) {
		u[i] = std::sqrt(x[i]);
	}
	evaluate_wave_function(u.data(), result, n, energy, lambda);
}

static const int SCksiODEdim = 3;

struct SCksiODEParams {
	gsl_spline   *density;
	gsl_interp_accel *acc;
	double e;
	double l;
};

int SCksiRHS(
	double       x, 
	const double y[], 
	double       dydx[], 
	void        *params
) {
	auto par = (SCksiODEParams *) params;

	auto  density = par->density;
	auto  acc = par->acc;
	auto& e	=   par->e;
	auto& l	=   par->l;

	double x2 = x*x;
	double p = std::abs(e*x2 + y[0] - l/x2);

	dydx[0] = 2.0*x*y[1];
	dydx[1] = x > 0 ? 2.0*gsl_spline_eval(density, x, acc)/x : 0.0;
	dydx[2] = -std::sqrt(p);

	return GSL_SUCCESS;
}

struct SCwfParams {
	gsl_spline        *wf;
	gsl_interp_accel *acc;
};

double sc_wf_norm(double x, void * params) {
	auto& p = *((SCwfParams*) params);
	double wf = gsl_spline_eval(p.wf, x, p.acc);
	return x*wf*wf;
}

void SemiclassicAtom::normalize_wave_function(
	const double* u, 
	double*       result, 
	std::size_t   n
) {
	gsl_spline *wfSpline;
	wfSpline = gsl_spline_alloc(gsl_interp_cspline, n);
	gsl_spline_init(wfSpline, u, result, n);
	auto wfacc = gsl_interp_accel_alloc();

	SCwfParams wfparams{wfSpline, wfacc};
	double wfNorm, wfError;
	gsl_function gslWF;
	gslWF.function = &sc_wf_norm;
	gslWF.params = &wfparams;

	gsl_integration_workspace *wrkSpace = gsl_integration_workspace_alloc(1024);
	gsl_integration_qag(
		&gslWF,            // function
		u[0],              // a
		u[n - 1],          // b
		tolerance,         // epsabs //? 0 or tolerance ?//
		tolerance,         // epsrel
		1024,              // limit
		GSL_INTEG_GAUSS41, // key
		wrkSpace,          // workspace
		&wfNorm,           // result
		&wfError           // abserr
	);

	wfNorm = std::sqrt(1.0/(2.0*r0*wfNorm));

	gsl_integration_workspace_free(wrkSpace);
	gsl_spline_free(wfSpline);
	gsl_interp_accel_free(wfacc);

	for (int i = 0; i < n; ++i) {
		result[i] *= wfNorm;
	}
}

void SemiclassicAtom::evaluate_wave_function(
	const double* u, 
	double*       result, 
	std::size_t   n,
	double        energy, 
	double        lambda
) {

	auto rpo = outerRP(energy, lambda);
	auto rpi = innerRP(energy, lambda);

	if (std::isnan(rpi) || std::isnan(rpo)) return;

	double umax = std::sqrt(rpo);
	double umin = std::sqrt(rpi);

	if (umax <= umin) return;

	double l = 0.5*lambda*lambda / r0 / r0;
	double e = energy;

	SCksiODEParams ksi_params;
	ksi_params.density = densSpline;
	ksi_params.acc = densAcc;
	ksi_params.e = e;
	ksi_params.l = l;

	double ksi0 = 0.0;

	double yksi[SCksiODEdim] = { 0 }; 

	gsl_odeiv2_step*	stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, SCksiODEdim);
	gsl_odeiv2_control* control = gsl_odeiv2_control_y_new (0.1*tolerance, 0.0);
	gsl_odeiv2_evolve*   evolve = gsl_odeiv2_evolve_alloc (SCksiODEdim);

	gsl_odeiv2_system sys = {SCksiRHS, NULL, SCksiODEdim, &ksi_params};

	if (umax < 1.0) {

		double from = 1.0;
		double to = umax;
		double h = -std::min(1.e-10, tolerance);

		while(from > to) {
			int status = gsl_odeiv2_evolve_apply (
				evolve, control, stepper,
				&sys,
				&from, to,
				&h, yksi
			);
			if (status != GSL_SUCCESS) {
				std::cerr << "ODE for action integration problems" << std::endl;
			}
		}

		ksi0 = 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);

		gsl_odeiv2_evolve_reset (evolve);
		gsl_odeiv2_step_reset (stepper);
	}

	double from = umax;
	double to = umin;
	yksi[2] = 0.0;
	double h = -tolerance;

	while(from > to) {
		int status = gsl_odeiv2_evolve_apply (
			evolve, control, stepper,
			&sys,
			&from, to,
			&h, yksi
		);
		if (status != GSL_SUCCESS) {
			std::cerr << "ODE for action integration problems" << std::endl;
		}
	}

	::aatk::util::specfunc::bessel::Jnu Jnu;
	::aatk::util::specfunc::bessel::Ynu Ynu;
	::aatk::util::specfunc::bessel::Knu Knu;

	double ksi21 = 1e-10 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);
	double Jp13 = Jnu(1.0/3.0, ksi21);
	double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21));
	double sign_2 = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;

	gsl_odeiv2_evolve_reset (evolve);
	gsl_odeiv2_step_reset (stepper);

	int i = n - 1;
	from = 1.0;
	yksi[0] = yksi[1] = yksi[2] = 0.0;
	h = -tolerance;

	std::vector<int> RPOinterpIdx;
	std::vector<int> RPIinterpIdx;
	RPOinterpIdx.reserve(256);
	RPIinterpIdx.reserve(256);

	double ksiRPOeps = 0.01;
	double ksiRPIeps = 0.1;

	while (u[i] > umax && i >= 0) {
		to = u[i];
		while(from > to) {
			int status = gsl_odeiv2_evolve_apply (
				evolve, control, stepper,
				&sys,
				&from, to,
				&h, yksi
			);
			if (status != GSL_SUCCESS) {
				std::cerr << "ODE for action integration problems" << std::endl;
			}
		}

		double ksi = 1e-10 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*yksi[2]);
		double p = std::sqrt(2.0*std::abs(e*to*to + yksi[0] - l / (to*to)))/to;

		if	    (ksi > 100.0) result[i] = 0.0;
		else if (ksi < ksiRPOeps) {
			RPOinterpIdx.push_back(i);
		}
		else {
			result[i] = sign_2/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);
		}
		--i; from = to;
	}

	gsl_odeiv2_evolve_reset (evolve);
	gsl_odeiv2_step_reset (stepper);

	from = umax;
	yksi[2] = 0.0;
	h = -tolerance;

	while (u[i] > umin && i >= 0) {

		to = u[i];
		while(from > to) {
			int status = gsl_odeiv2_evolve_apply (
				evolve, control, stepper,
				&sys,
				&from, to,
				&h, yksi
			);
			if (status != GSL_SUCCESS) {
				std::cerr << "ODE for action integration problems" << std::endl;
			}
		}

		double ksi2x  = 1e-10 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]); 
		double a = ksi2x/ksi21;

		double p = std::sqrt(2.0*std::abs(e*to*to + yksi[0] - l / (to*to)))/to;

		double ksi1x  = 1e-10 + std::abs(ksi21 - ksi2x);

		double Jp13_1 = Jnu(1.0/3.0, ksi1x);
		double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

		double Jp13_2 = Jnu(1.0/3.0, ksi2x);
		double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

		double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
		double R2 = sign_2*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

		result[i] = a*R1 + (1.0 - a)*R2;

		if (umax > 1.0 - tolerance) {
			result[i] = R1;
		}
		if (ksi2x < ksiRPOeps && umax <= 1.0 - tolerance) {
			RPOinterpIdx.push_back(i);
		}
		if (ksi1x < ksiRPIeps) {
			RPIinterpIdx.push_back(i);
		}

		--i; from = to;
	}

	gsl_odeiv2_evolve_reset (evolve);
	gsl_odeiv2_step_reset (stepper);

	from = umin;
	yksi[2] = 0.0;
	h = -std::min(1.e-10, tolerance);

	while (u[i] < umin && i >= 0) {

		double to = u[i];

		if (u[i]*u[i] > tolerance) {

			while(from > to) {
				int status = gsl_odeiv2_evolve_apply (
					evolve, control, stepper,
					&sys,
					&from, to,
					&h, yksi
				);
				if (status != GSL_SUCCESS) {
					std::cerr << "ODE for action integration problems" << std::endl;
				}
			}

			double ksi = 1e-10 + 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);
			double p = std::sqrt(2.0*std::abs(e*to*to + yksi[0] - l / (to*to)))/to;

			if	    (ksi > 100.0) result[i] = 0.0;
			else if (ksi < ksiRPIeps) {
				RPIinterpIdx.push_back(i);
			}
			else {
				result[i] = 1.0/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);
			}
		}
		else { result[i] = 0.0; }

		--i; from = to;
	}

	// interpolation around return points
	if (RPOinterpIdx.size() > 0) {
		int iR = RPOinterpIdx.front() + 1;
		int iL = RPOinterpIdx.back() - 1;
		double wfR = result[iR];
		double wfL = result[iL];
		double xL = std::pow(u[iL], 4.0/3.0);
		double xR = std::pow(u[iR], 4.0/3.0);
		for (int i = iR - 1; i > iL; --i) {
			double x = std::pow(u[i], 4.0/3.0);
			result[i] = wfL + (wfR - wfL)/(xR - xL)*(x - xL);
		}
	}
	if (RPIinterpIdx.size() > 0) {
		int iR = RPIinterpIdx.front() + 1;
		int iL = RPIinterpIdx.back() - 1;
		double wfR = result[iR];
		double wfL = result[iL];
		double xL = std::pow(u[iL], 4.0/3.0);
		double xR = std::pow(u[iR], 4.0/3.0);
		for (int i = iR - 1; i > iL; --i) {
			double x = std::pow(u[i], 4.0/3.0);
			result[i] = wfL + (wfR - wfL)/(xR - xL)*(x - xL);
		}
	}

	normalize_wave_function(u, result, n);
}

}
}