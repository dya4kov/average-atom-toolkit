#include <cmath>
#include <algorithm>
#include <iostream>
// #include <functional>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

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

static const int SCwaveFunctionNormODEdim = 4;

struct SCwaveFunctionNormODEParams {
	gsl_spline   *density;
	gsl_interp_accel *acc;
	double e;
	double l;
	double r0;
	double ksi0;
	double ksi21;
	double xmin;
	double xmax;
	double sign_2;
	aatk::util::specfunc::bessel::Jnu Jnu;
	aatk::util::specfunc::bessel::Ynu Ynu;
	aatk::util::specfunc::bessel::Knu Knu;
};

int SCwaveFunctionNormRHS(
	double       x, 
	const double y[], 
	double       dydx[], 
	void        *params
) {

	auto par = (SCwaveFunctionNormODEParams *) params;

	auto  density = par->density;
	auto  acc     = par->acc;

	auto& e	      = par->e;
	auto& l	      = par->l;
	auto& r0      = par->r0;
	auto& ksi0    = par->ksi0;
	auto& ksi21   = par->ksi21;
	auto& xmin    = par->xmin;
	auto& xmax    = par->xmax;
	auto& sign_2  = par->sign_2;

	auto& Jnu     = par->Jnu;
	auto& Ynu     = par->Ynu;
	auto& Knu     = par->Knu;

	double p2half = std::abs(e*x*x + y[0] - l/(x*x));
	double p      = std::sqrt(2.0*p2half);

	dydx[0] = 2.0*x*y[1];
	dydx[1] = x > 0 ? 2.0*gsl_spline_eval(density, x, acc)/x : 0.0;
	dydx[2] = -std::sqrt(p2half);

	if (x > xmax) {
		double ksi = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2]));
		// std::cout << "x = " << x << ", xmax = " << xmax << ", ksi0 = " << ksi0 << ", ksi = " << ksi << ", p = " << p << ", y[3] = " << y[3] << std::endl;
		dydx[3] =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
		dydx[3] = -x*x*ksi/p*dydx[3]*dydx[3];
	}

	if (x > xmin && x <= xmax) {

		double ksi2x = 1e-8 + std::abs(ksi0 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2])); 
		double ksi1x = 1e-8 + std::abs(ksi21 - ksi2x);
		double ax   = ksi2x/ksi21;

		// std::cout << "x = " << x << ", ksi21 = " << ksi21 << ", ksi2x = " << ksi2x << ", ksi1x = " << ksi1x << ", p = " << p << ", y[3] = " << y[3] << std::endl;
		
		double Jp13_1 = Jnu(1.0/3.0, ksi1x);
		double Jm13_1 = 0.5*(Jp13_1 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi1x));

		double Jp13_2 = Jnu(1.0/3.0, ksi2x);
		double Jm13_2 = 0.5*(Jp13_2 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi2x));

		double R1 = std::sqrt(ksi1x/(3.0*p))*(Jp13_1 + Jm13_1);
		double R2 = sign_2*std::sqrt(ksi2x/(3.0*p))*(Jp13_2 + Jm13_2);

		dydx[3] = xmax > 1.0 - 1e-8 ? R1 : (ax*R1 + (1.0 - ax)*R2);
		dydx[3] *= -x*x*dydx[3];

		// std::cout << "x = " << x << " Jp13_2 = " << Jp13_2 << " Jm13_2 = " << Jm13_2 << " ksi2x = " << ksi2x << " Ynu(1.0/3.0, ksi2x) " << gsl_sf_bessel_Ynu(1.0/3.0, ksi2x) << ", R1 = " << R1 << ", R2 = " << R2 << ", R1 + R2 = " << ax*R1 + (1.0 - ax)*R2 << ", dydx[3] = " << dydx[3] << std::endl;
	}

	if (x <= xmin) {
		double ksi = 1e-8 + std::abs(ksi0 + ksi21 - 2.0*r0*std::sqrt(2.0)*std::abs(y[2]));

		// std::cout << "x = " << x << ", ksi = " << ksi << ", p = " << p << ", y[3] = " << y[3] << std::endl;

		dydx[3] =  ksi > 100.0 ? 0.0 : Knu(1.0/3.0, ksi)/M_PI;
		dydx[3] = -x*x*ksi/p*dydx[3]*dydx[3];
	}

	return GSL_SUCCESS;
}

double SemiclassicAtom::waveFunctionNorm(
	double energy, 
	double lambda, 
	double rpi, 
	double rpo
) {

	double umax = std::sqrt(rpo);
	double umin = std::sqrt(rpi);

	if (umax <= umin) return 1.0;

	double l = 0.5*lambda*lambda / r0 / r0;
	double e = energy;

	SCksiODEParams ksi_params;
	ksi_params.density = densSpline;
	ksi_params.acc = acc;
	ksi_params.e = e;
	ksi_params.l = l;

	double yksi[SCksiODEdim] = { 0 }; 

	gsl_odeiv2_step*	stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, SCksiODEdim);
	gsl_odeiv2_control* control = gsl_odeiv2_control_y_new (0.1*tolerance, 0.0);
	gsl_odeiv2_evolve*   evolve = gsl_odeiv2_evolve_alloc (SCksiODEdim);

	gsl_odeiv2_system sys = {SCksiRHS, NULL, SCksiODEdim, &ksi_params};

	double ksi0 = 0.0;

	if (umax < 1.0) {

		double from = 1.0;
		double to = umax;
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

		ksi0 = 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);

		gsl_odeiv2_evolve_reset (evolve);
		gsl_odeiv2_step_reset (stepper);
	}

	double from = umax;
	double to = umin;
	// yksi[0] = gsl_spline_eval(phiSpline, from, acc);
	// yksi[1] = gsl_spline_eval(dphiSpline, from, acc);
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

	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (stepper);

	::aatk::util::specfunc::bessel::Jnu Jnu;
	::aatk::util::specfunc::bessel::Ynu Ynu;
	::aatk::util::specfunc::bessel::Knu Knu;

	double ksi21 = 2.0*r0*std::sqrt(2.0)*std::abs(yksi[2]);
	double Jp13 = Jnu(1.0/3.0, ksi21);
	double Jm13 = 0.5*(Jp13 - std::sqrt(3.0)*Ynu(1.0/3.0, ksi21));
	double sign_2 = Jp13 + Jm13 > 0.0 ? 1.0 : -1.0;

	SCwaveFunctionNormODEParams wf_params;
	wf_params.e       = e;
	wf_params.l       = l;
	wf_params.r0      = r0;
	wf_params.density = densSpline;
	wf_params.acc     = acc;
	wf_params.ksi0    = ksi0;
	wf_params.ksi21   = ksi21;
	wf_params.xmin    = umin;
	wf_params.xmax    = umax;
	wf_params.sign_2  = sign_2;

	stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, SCwaveFunctionNormODEdim);
	control = gsl_odeiv2_control_y_new (0.1*tolerance, 0.0);
	evolve  = gsl_odeiv2_evolve_alloc (SCwaveFunctionNormODEdim);

	sys = {SCwaveFunctionNormRHS, NULL, SCwaveFunctionNormODEdim, &wf_params};

	from = 1.0;
	to = std::sqrt(0.1*tolerance);
	double ynorm[SCwaveFunctionNormODEdim] = { 0 }; 
	h = -tolerance;

	while(from > to) {
		int status = gsl_odeiv2_evolve_apply (
			evolve, control, stepper,
			&sys,
			&from, to,
			&h, ynorm
		);
		if (status != GSL_SUCCESS) {
			std::cerr << "ODE for wave function norm integration problems" << std::endl;
		}
	}

	return std::sqrt(1.0/(2.0*r0*std::abs(ynorm[3])));
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

	double normValue = waveFunctionNorm(energy, lambda, rpi, rpo);

	double umax = std::sqrt(rpo);
	double umin = std::sqrt(rpi);

	if (umax <= umin) return;

	double l = 0.5*lambda*lambda / r0 / r0;
	double e = energy;

	SCksiODEParams ksi_params;
	ksi_params.density = densSpline;
	ksi_params.acc = acc;
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
	// yksi[0] = gsl_spline_eval(phiSpline, from, acc);
	// yksi[1] = gsl_spline_eval(dphiSpline, from, acc);
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

	// gamma(1/3)
	double gamma   = 2.6789385347077476336557;

	double RPOx = umax*umax;
	double RPOphi = gsl_spline_eval(phiSpline, umax, acc);
	double RPOdphi = gsl_spline_eval(dphiSpline, umax, acc);
	double RPOc = 2.0/RPOx*((RPOdphi - RPOphi/RPOx) + 2.0*l/(RPOx*RPOx));

	double ksi16pm12_2 = std::pow(2.0*r0/(3.0*std::abs(RPOc)), 1.0/6.0);

	double R2approx = sign_2*normValue*gamma/M_PI*ksi16pm12_2*std::pow(2.0, -2.0/3.0);

	double RPIx = umin*umin;
	double RPIphi = gsl_spline_eval(phiSpline, umin, acc);
	double RPIdphi = gsl_spline_eval(dphiSpline, umin, acc);
	double RPIc = 2.0/RPIx*((RPIdphi - RPIphi/RPIx) + 2.0*l/(RPIx*RPIx));

	double ksi16pm12_1 = std::pow(2.0*r0/(3.0*std::abs(RPIc)), 1.0/6.0);

	double R1approx = normValue*gamma/M_PI*ksi16pm12_1*std::pow(2.0, -2.0/3.0);

	double ksiParam = std::pow(2.0, 1.0/3.0)*std::sqrt(3.0)*M_PI/(gamma*gamma);
	double ksiRPeps = 2e-3;

	gsl_odeiv2_evolve_reset (evolve);
	gsl_odeiv2_step_reset (stepper);

	int i = n - 1;
	from = 1.0;
	yksi[0] = yksi[1] = yksi[2] = 0.0;
	h = -tolerance;

	std::vector<int> RPOlinearInterpIdx;
	std::vector<int> RPIlinearInterpIdx;
	RPOlinearInterpIdx.reserve(50);
	RPIlinearInterpIdx.reserve(50);

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
		else if (ksi < ksiRPeps) {
			// result[i] = R2approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
			RPOlinearInterpIdx.push_back(i);
		}
		else {
			result[i] = sign_2*normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);
		}

		--i; from = to;
	}

	gsl_odeiv2_evolve_reset (evolve);
	gsl_odeiv2_step_reset (stepper);

	from = umax;
	// yksi[0] = gsl_spline_eval(phiSpline, from, acc);
	// yksi[1] = gsl_spline_eval(dphiSpline, from, acc);
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

		result[i] = normValue*(a*R1 + (1.0 - a)*R2);

		if (umax > 1.0 - tolerance) {
			result[i] = normValue*R1;
		}
		if (ksi2x < ksiRPeps && umax <= 1.0 - tolerance) {
			RPOlinearInterpIdx.push_back(i);
		}
		if (ksi1x < ksiRPeps) {
			RPIlinearInterpIdx.push_back(i);
		}

		// if (ksi2x < ksiRPeps && umax < 1.0 - tolerance) {
		// 	result[i] = R2approx*(1.0 + ksiParam*std::pow(ksi2x, 2.0/3.0));
		// }
		// if (ksi1x < ksiRPeps) {
		// 	result[i] = R1approx*(1.0 + ksiParam*std::pow(ksi1x, 2.0/3.0));
		// }

		--i; from = to;
	}

	gsl_odeiv2_evolve_reset (evolve);
	gsl_odeiv2_step_reset (stepper);

	from = umin;
	// yksi[0] = gsl_spline_eval(phiSpline, from, acc);
	// yksi[1] = gsl_spline_eval(dphiSpline, from, acc);
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
			else if (ksi < ksiRPeps) {
				// result[i] = R1approx*(1.0 - ksiParam*std::pow(ksi, 2.0/3.0));
				RPIlinearInterpIdx.push_back(i);
			}
			else {
				result[i] = normValue/M_PI*std::sqrt(ksi/p)*Knu(1.0/3.0, ksi);
			}
		}
		else { result[i] = 0.0; }

		--i; from = to;
	}

	// linear interpolation around return points
	if (RPOlinearInterpIdx.size() > 0) {
		int iR = RPOlinearInterpIdx.front() + 1;
		int iL = RPOlinearInterpIdx.back() - 1;
		double wfR = result[iR];
		double wfL = result[iL];
		double uR = u[iR];
		double uL = u[iL];
		for (int i = iR - 1; i > iL; --i) {
			result[i] = wfR - (wfR - wfL)/(uR - uL)*(uR - u[i]);
			std::cout << i << "  " << u[i]*u[i] << " " << result[i] << std::endl;
		}
	}
	if (RPIlinearInterpIdx.size() > 0) {
		int iR = RPIlinearInterpIdx.front() + 1;
		int iL = RPIlinearInterpIdx.back() - 1;
		double wfR = result[iR];
		double wfL = result[iL];
		double uR = u[iR];
		double uL = u[iL];
		for (int i = iR - 1; i > iL; --i) {
			result[i] = wfR - (wfR - wfL)/(uR - uL)*(uR - u[i]);
			std::cout << i << "  " << u[i]*u[i] << " " << result[i] << std::endl;
		}
	}
}

}
}