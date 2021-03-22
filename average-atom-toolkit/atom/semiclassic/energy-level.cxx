#include <cmath>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

#include <average-atom-toolkit/atom/semiclassic.h>

namespace aatk {
namespace atom {

struct SCRPparams {
	gsl_spline       *phi;
    gsl_interp_accel *acc;
    double e;
    double l;
};

double SCRPfunc(double x, void *params) {
	auto p = (SCRPparams *) params;
	double phi = gsl_spline_eval(p->phi, x, p->acc);
	double p2 = p->e*x*x + phi - p->l / (x*x);
	return -p2;
}

double SemiclassicAtom::innerRP(double energy, double lambda) {
	double l = 0.5*lambda*lambda / r0 / r0;
	double e = energy;

	SCRPparams params;
	params.phi = phiSpline;
	params.acc = phiAcc;
	params.e = e;
	params.l = l;

	double uLeft = 1.e-6;
	double uMid = 0.9;
	double uRight = 1.0;

	gsl_function F;
	F.function = &SCRPfunc;
	F.params = &params;

	int status;
	int iter, max_iter = 100;

	if (e < l) {
		// check that maximum exists
		double p2right = -SCRPfunc(uRight, &params);
		double p2mid = -SCRPfunc(uMid, &params);
		if (p2mid > p2right) {
			// locate maximum of semiclassic momentum
			gsl_min_fminimizer* minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
			gsl_min_fminimizer_set(minimizer, &F, uMid, uLeft, uRight);
			iter = 0;
			status = GSL_CONTINUE;
			while (status == GSL_CONTINUE && iter < max_iter) {
				status = gsl_min_fminimizer_iterate(minimizer);
				uLeft  = gsl_min_fminimizer_x_lower(minimizer);
				uRight = gsl_min_fminimizer_x_upper(minimizer);
				status = gsl_min_test_interval(uLeft, uRight, 0.0, 0.1*tolerance);
			}
			gsl_min_fminimizer_free(minimizer);
			uRight = uLeft;
		}
	}
	// check that root exists
	double p2 = -SCRPfunc(uRight, &params);
	if (p2 < 0.0) return NAN; // no return point
	// locate root of semiclassic momentum
	uLeft = 1.e-6;
	gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	gsl_root_fsolver_set(solver, &F, uLeft, uRight);
	iter = 0;
	status = GSL_CONTINUE;
	while (status == GSL_CONTINUE && iter < max_iter) {
		status = gsl_root_fsolver_iterate (solver);
		uLeft  = gsl_root_fsolver_x_lower (solver);
		uRight = gsl_root_fsolver_x_upper (solver);
		status = gsl_root_test_interval(uLeft, uRight, 0.0, 0.1*tolerance);
	}
	gsl_root_fsolver_free(solver);
	return uRight*uRight;
}

double SemiclassicAtom::outerRP(double energy, double lambda) {

	double l = 0.5*lambda*lambda / r0 / r0;
	double e = energy;

	if (e >= l) return 1.0;

	SCRPparams params;
	params.phi = phiSpline;
	params.acc = phiAcc;
	params.e = e;
	params.l = l;

	double uLeft = 1.e-6;
	double uMid = 0.9;
	double uRight = 1.0;

	gsl_function F;
	F.function = &SCRPfunc;
	F.params = &params;

	int status;
	int iter, max_iter = 100;

	// check that maximum exists
	double p2right = -SCRPfunc(uRight, &params);
	double p2mid = -SCRPfunc(uMid, &params);
	if (p2mid > p2right) {
		// locate maximum of semiclassic momentum
		gsl_min_fminimizer* minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
		gsl_min_fminimizer_set(minimizer, &F, uMid, uLeft, uRight);
		iter = 0;
		status = GSL_CONTINUE;
		while (status == GSL_CONTINUE && iter < max_iter) {
			status = gsl_min_fminimizer_iterate(minimizer);
			uLeft  = gsl_min_fminimizer_x_lower(minimizer);
			uRight = gsl_min_fminimizer_x_upper(minimizer);
			status = gsl_min_test_interval(uLeft, uRight, 0.0, 0.1*tolerance);
		}
		gsl_min_fminimizer_free(minimizer);
		uRight = uLeft;
	}
	uLeft = uRight;
	uRight = 1.0;
	// check that root exists
	double p2 = -SCRPfunc(uLeft, &params);
	if (p2 < 0.0) return NAN; // no return point
	// locate root of semiclassic momentum
	gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	gsl_root_fsolver_set(solver, &F, uLeft, uRight);
	iter = 0;
	status = GSL_CONTINUE;
	while (status == GSL_CONTINUE && iter < max_iter) {
		status = gsl_root_fsolver_iterate (solver);
		uLeft  = gsl_root_fsolver_x_lower (solver);
		uRight = gsl_root_fsolver_x_upper (solver);
		status = gsl_root_test_interval(uLeft, uRight, 0.0, 0.1*tolerance);
	}
	gsl_root_fsolver_free(solver);
	return uLeft*uLeft;
}

static const int SCactionODEdim = 3;

struct SCactionODEParams {
    gsl_spline   *density;
    gsl_interp_accel *acc;
    double e;
    double l;
};

int SCactionRHS(double x, const double y[], double dydx[], void *params) {
    auto p = (SCactionODEParams *) params;

    dydx[0] = 2.0*x*y[1];
    dydx[1] = x > 0 ? 2.0/x*gsl_spline_eval(p->density, x, p->acc) : 0.0;
    double p2 = p->e*x*x + y[0] - p->l/(x*x);
    dydx[2] = p2 > 0.0 ? -std::sqrt(p2) : 0.0;

    return GSL_SUCCESS;
}

double SemiclassicAtom::action(double energy, double lambda) {
	auto rpi = innerRP(energy, lambda);
	auto rpo = outerRP(energy, lambda);

	if (std::isnan(rpi) || std::isnan(rpo)) return 0.0;

    double umin = std::sqrt(rpi);
	double umax = std::sqrt(rpo);

    double result = 0.0;

    if (umax <= umin) return result;

	double l = 0.5*lambda*lambda / r0 / r0;
	double e = energy;

    SCactionODEParams params;
	params.density = densSpline;
	params.acc = densAcc;
	params.e = e;
	params.l = l;

	double ay[SCactionODEdim] = { 0 }; 

	gsl_odeiv2_step*    stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, SCactionODEdim);
	gsl_odeiv2_control* control = gsl_odeiv2_control_y_new (0.1*tolerance, 0.0);
	gsl_odeiv2_evolve*   evolve = gsl_odeiv2_evolve_alloc (SCactionODEdim);

	gsl_odeiv2_system sys = {SCactionRHS, NULL, SCactionODEdim, &params};

	double from = umax;
	double to = umin;
	ay[0] = gsl_spline_eval(phiSpline, umax, phiAcc);
	ay[1] = gsl_spline_eval(dphiSpline, umax, dphiAcc);
	double h = -std::min(1.e-10, tolerance);

	while(from > to) {
		int status = gsl_odeiv2_evolve_apply (
			evolve, control, stepper,
			&sys,
			&from, to,
			&h, ay
		);
		if (status != GSL_SUCCESS) {
			std::cerr << "ODE for action integration problems" << std::endl;
		}
	}

	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (stepper);

	return 2.0*r0*std::sqrt(2.0)*ay[2];
}

double SemiclassicAtom::energyLevel(int n, int l) {
	if (!eLevelReady[n][l]) evaluate_energy_level(n, l);
	return eLevel[n][l];
}

struct SCeLevelParams {
	SemiclassicAtom *atom;
    double exactAction;
    double lambda;
};

double SCeLevelfunc(double e, void *params) {
	auto p = (SCeLevelParams *) params;
	double result = p->atom->action(e, p->lambda) - p->exactAction;
	return result;
}

void SemiclassicAtom::evaluate_energy_level(int n, int l) {
    double exact = M_PI*(n - l - 0.5);
    double lambda = l + 0.5;

    int ie = 0;
    double act = -1;
    while (act < exact && ie < eLevelStart.size()) {
    	act = action(eLevelStart[ie], lambda); ++ie;
    }

    SCeLevelParams params;
	params.atom = this;
	params.exactAction = exact;
	params.lambda = lambda;

    double eMin = eLevelStart[ie - 2];
    double eMax = eLevelStart[ie - 1];

	gsl_function F;
	F.function = &SCeLevelfunc;
	F.params = &params;

	int status;
	int iter, max_iter = 100;

	// locate root for Bohr-Sommerfeld condition
	gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	gsl_root_fsolver_set(solver, &F, eMin, eMax);
	iter = 0;
	status = GSL_CONTINUE;
	while (status == GSL_CONTINUE && iter < max_iter) {
		status = gsl_root_fsolver_iterate (solver);
		eMin   = gsl_root_fsolver_x_lower (solver);
		eMax   = gsl_root_fsolver_x_upper (solver);
		status = gsl_root_test_interval(eMin, eMax, 0.0, tolerance);
	}
	gsl_root_fsolver_free(solver);
    
    eLevel[n][l] = 0.5*(eMin + eMax);
    eLevelReady[n][l] = true;
}

void SemiclassicAtom::evaluate_boundary_energy(){
    int n, n_border;
    double E_curr = 0;
    bool check = false;

    for ( n = 1; n <= Z && !check ; ++n) {
        if (n > nmax){
            nmax++;
            std::vector<double> l_vector(nmax);
            std::vector<bool> l_bool_vector(nmax);
            eLevel.push_back(l_vector);
            eLevelReady.push_back(l_bool_vector);
        }

        for (int l = 0; l < n; ++l) {
            evaluate_energy_level(n, l);
            E_curr = eLevel[n][l];
            if (E_curr >= 0 && !check){
                boundaryEnergy = E_curr;
                check = true;
                n_border = n;
            }
        }
    }

    check = false;

    for ( n = n_border + 1; n <= Z && !check ; ++n) {
        if (n > nmax){
            nmax++;
            std::vector<double> l_vector(nmax);
            std::vector<bool> l_bool_vector(nmax);
            eLevel.push_back(l_vector);
            eLevelReady.push_back(l_bool_vector);
        }

        check = true;

        for (int l = 0; l < n; ++l) {
            evaluate_energy_level(n, l);
            E_curr = eLevel[n][l];
            if (E_curr >=0 && E_curr < boundaryEnergy) {
                boundaryEnergy = E_curr;
            }
            check = check && (E_curr > 0);
        }
    }

    nmax = n - 2 ; // ?check 
}

}
}