#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

#include <average-atom-toolkit/transport/Kfunction.h>

namespace aatk {
namespace transport {

static const double max_exponent = 100.;

struct K_EqParams {
    double  n;
    double  T;
    double  M;
    Spline* tau;
};

int K_odeFunction(double e, const double y[], double f[], void *params) {
    K_EqParams* p = (K_EqParams *)params;
    double n   = p->n;
    double T   = p->T;
    double M   = p->M;
    auto&  tau = *(p->tau);
    double exponent = (e - M)/T;
    if (std::abs(exponent) < max_exponent)
    	 f[0] = std::pow(e, n)*tau(std::log10(e))/(1. + std::cosh(exponent));
    else f[0] = 0.;
    return GSL_SUCCESS;
}

Kfunction::Kfunction(double _n, double _T, double _M) :
                          n(_n),     T(_T),     M(_M)
{}

Kfunction::~Kfunction() {
	delete interpolation;
}

double Kfunction::operator()(double* _e, double* _tau, std::size_t npoints) {
	std::vector<double> e(npoints);
    std::vector<double> tau(npoints);
	bool needSort = false;
	for (std::size_t i = 0; i < npoints - 1; ++i) {
		needSort = needSort || _e[i] > _e[i + 1];
	}
	if (needSort) {
		std::vector<std::size_t> idx(npoints);
    	std::iota(idx.begin(), idx.end(), 0);
    	std::sort(idx.begin(), idx.end(),
    	   [&_e](std::size_t i1, std::size_t i2) {
    	        return _e[i1] < _e[i2];
    	   }
    	);
    	std::size_t k = 0;
    	for (auto i : idx) {
			e[k] = std::log10(_e[i]);
			tau[k] = _tau[i];
			++k;
    	}
	}
	else {
		for (std::size_t i = 0; i < npoints; ++i) {
			e[i] = std::log10(_e[i]);
			tau[i] = _tau[i];
    	}
	}
    interpolation = new Spline(e, tau);
	auto& tau_func = *interpolation;

	emin = std::pow(10., e.front());
	emax = std::pow(10., e.back());

	const int dim = 1;
    K_EqParams params;
    params.n = n;
    params.T = T;
    params.M = M;
    params.tau = interpolation;

    const gsl_odeiv2_step_type * stepT
            = gsl_odeiv2_step_rk8pd;

    gsl_odeiv2_step * stepper
            = gsl_odeiv2_step_alloc (stepT, dim);
    gsl_odeiv2_control * control
            = gsl_odeiv2_control_y_new (0.0, 1.e-6);
    gsl_odeiv2_evolve * evolve
            = gsl_odeiv2_evolve_alloc (dim);

    gsl_odeiv2_system sys = {K_odeFunction, NULL, dim, &params};

    double from   = emin;
	double to     = emax;
    double h      = 1.e-6;
    double y[dim] = { 0 };

    while(from < to) {
        int status = gsl_odeiv2_evolve_apply (evolve, control, stepper,
                                              &sys,
                                              &from, to,
                                              &h, y);
        // std::cout << "e = " << from << ", h = " << h << ", K = " << y[0] << ", (e - M)/T = " << (from - M)/T << std::endl;
        if (status != GSL_SUCCESS) {
            std::cout << "ODE for Kfunction integration problems" << std::endl;
        }
    }
    gsl_odeiv2_evolve_free (evolve);
    gsl_odeiv2_control_free (control);
    gsl_odeiv2_step_free (stepper);
    return -0.5*y[0]/T;
}

}
}
