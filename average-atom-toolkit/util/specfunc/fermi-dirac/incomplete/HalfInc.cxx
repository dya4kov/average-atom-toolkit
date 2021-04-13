#include <cmath>
#include <algorithm>
#include <iostream>
#include <average-atom-toolkit/util/specfunc/fermi-dirac/incomplete.h>

namespace aatk {
namespace util {
namespace specfunc {
namespace FDI {

static const int FDHalfIncODEdim = 1;

int FDHalfIncRHS(double x, const double y[], double dydx[], void *params) {
	auto& p = *((double *) params);
	dydx[0] = std::sqrt(x)/(1.0 + std::exp(x - p));
	return GSL_SUCCESS;
}

HalfInc::HalfInc(double tolerance) {
	stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, FDHalfIncODEdim);
	control = gsl_odeiv2_control_y_new (0.0, tolerance);
	evolve = gsl_odeiv2_evolve_alloc (FDHalfIncODEdim);
}

HalfInc::~HalfInc() {
	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (stepper);
}

double HalfInc::value(const double x, const double y) {
	double offset = 25.0;
	double result = 0.0;

	if (x - offset > 0.0) {
		result += 2.0/3.0*std::min(y, x - offset)*std::sqrt(std::min(y, x - offset));
	}
	if (y > x - offset) {
		double yode[FDHalfIncODEdim] = { 0 }; 
		double params = x;
		gsl_odeiv2_system sys = {FDHalfIncRHS, NULL, FDHalfIncODEdim, &params};
		double from = std::max(0.0, x - offset);
		double to = std::min(y, x + offset);
		double h = std::min(1.e-6, tolerance);
		while(from < to) {
			int status = gsl_odeiv2_evolve_apply (
				evolve, control, stepper,
				&sys,
				&from, to,
				&h, yode
			);
			if (status != GSL_SUCCESS) {
				std::cerr << "ODE for FDHalfInc integration problems" << std::endl;
			}
		}
		result += yode[0];
		gsl_odeiv2_evolve_reset(evolve);
		gsl_odeiv2_step_reset(stepper);
	}
	if (y > x + offset) {
		result += std::exp(x)*(
		  0.5*std::sqrt(M_PI)*(
			  std::erf(std::sqrt(y)) - 
			  std::erf(std::sqrt(std::max(0.0, x + offset)))
		  )
		  -
		  (
			  std::sqrt(y)*std::exp(-y) -
			std::sqrt(std::max(0.0, x + offset))*std::exp(-std::max(0.0, x + offset))
		  )
		);
	}

	return result;
}

}
}
}
}