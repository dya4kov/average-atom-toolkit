#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

#include <average-atom-toolkit/atom/semiclassic.h>

namespace aatk {
namespace atom {

static const int SCpotentialODEdim = 2;

struct SCpotentialODEParams {
	gsl_spline   *density;
	gsl_interp_accel *acc;
};

int SCpotentialRHS(double x, const double y[], double dydx[], void *params) {
	auto p = (SCpotentialODEParams *) params;

	dydx[0] = 2.0*x*y[1];
	dydx[1] = x > 0 ? 2.0/x*gsl_spline_eval(p->density, x, p->acc) : 0.0;

	return GSL_SUCCESS;
}

void SemiclassicAtom::evaluate_potential() {
	// u - sorted array from 0.0 to 1.0
	const double* u = mesh.data();
	double*       y = pot.data();
	double*      dy = dpot.data();

	SCpotentialODEParams params;
	params.density = densSpline;
	params.acc = densAcc;

	double phi[2] = { 0 }; 

	gsl_odeiv2_step*    stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, SCpotentialODEdim);
	gsl_odeiv2_control* control = gsl_odeiv2_control_y_new (0.0, 0.1*tolerance);
	gsl_odeiv2_evolve*   evolve = gsl_odeiv2_evolve_alloc (SCpotentialODEdim);

	gsl_odeiv2_system sys = {SCpotentialRHS, NULL, SCpotentialODEdim, &params};

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

}
}