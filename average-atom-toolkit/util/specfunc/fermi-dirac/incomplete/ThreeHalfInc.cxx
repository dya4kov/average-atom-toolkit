#include <cmath>
#include <algorithm>
#include <iostream>
#include <average-atom-toolkit/util/specfunc/fermi-dirac/incomplete.h>

namespace aatk {
namespace util {
namespace specfunc {
namespace FDI {

static const int FDThreeHalfIncODEdim = 1;

int FDThreeHalfIncRHS(double x, const double y[], double dydx[], void *params) {
	auto& p = *((double *) params);
	dydx[0] = x*std::sqrt(x)/(1.0 + std::exp(x - p));
	return GSL_SUCCESS;
}

ThreeHalfInc::ThreeHalfInc(double tolerance) {
	stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, FDThreeHalfIncODEdim);
	control = gsl_odeiv2_control_y_new (0.0, 0.1*tolerance);
	evolve  = gsl_odeiv2_evolve_alloc (FDThreeHalfIncODEdim);
}

ThreeHalfInc::~ThreeHalfInc() {
	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (stepper);
}

double ThreeHalfInc::value(const double x, const double y) {
	double result = 0;
	double yode[FDThreeHalfIncODEdim] = { 0 }; 
	double params = x;
	gsl_odeiv2_system sys = {FDThreeHalfIncRHS, NULL, FDThreeHalfIncODEdim, &params};
	double from = 0.0;
	double to = y;
	double h = std::min(1.e-6, tolerance);
	while(from < to) {
		int status = gsl_odeiv2_evolve_apply (
			evolve, control, stepper,
			&sys,
			&from, to,
			&h, yode
		);
		if (status != GSL_SUCCESS) {
			std::cerr << "ODE for FDThreeHalfInc integration problems" << std::endl;
		}
	}
	result = yode[0];
	gsl_odeiv2_evolve_reset(evolve);
	gsl_odeiv2_step_reset(stepper);
    return result;
}

}
}
}
}