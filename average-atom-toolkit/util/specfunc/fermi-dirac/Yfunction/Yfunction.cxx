#include <cmath>
#include <algorithm>
#include <iostream>
#include <average-atom-toolkit/util/specfunc/fermi-dirac/Yfunction.h>

using namespace aatk::util::specfunc;

static const int YfunctionODEdim = 1;

struct YfunctionODEParams {
    FermiDirac<MHalf> FDmhalf;
};

int Yfunction_rhs(double x, const double y[], double dydx[], void *params) {
	auto& p = *((YfunctionODEParams *) params);
	auto& FDmhalf = p.FDmhalf;
	double f = FDmhalf(x);
	dydx[0] = f*f;
	return GSL_SUCCESS;
}

Yfunction::Yfunction(double _tolerance) : tolerance(_tolerance) {
	stepper = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rk8pd, YfunctionODEdim);
	control = gsl_odeiv2_control_y_new (0.0, 0.1*tolerance);
	evolve  = gsl_odeiv2_evolve_alloc (YfunctionODEdim);
}

Yfunction::~Yfunction() {
	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (stepper);
}

double Yfunction::operator()(const double x) {
	return 0.5*FDhalf(x)*FDmhalf(x) + 1.5*integral(x);
}

double Yfunction::derivative(const double x) {
	double f = FDmhalf(x);
	return 1.75*f*f + 0.5*FDhalf(x)*FDdmhalf(x);
}

std::vector<double>& Yfunction::operator()(const std::vector<double>& x) {
	size_t n = x.size();
	auto result = new std::vector<double>(n);
	for (size_t i = 0; i < n; ++i) {
		(*result)[i] = operator()(x[i]);
	}
	return *result;
}
std::vector<double>& Yfunction::derivative(const std::vector<double>& x) {
	size_t n = x.size();
	auto result = new std::vector<double>(n);
	for (size_t i = 0; i < n; ++i) {
		(*result)[i] = derivative(x[i]);
	}
	return *result;
}

double* Yfunction::operator()(const double* x, const size_t& n) {
	auto result = new double[n];
	for (size_t i = 0; i < n; ++i) {
		result[i] = operator()(x[i]);
	}
	return result;
}

double* Yfunction::derivative(const double* x, const size_t& n) {
	auto result = new double[n];
	for (size_t i = 0; i < n; ++i) {
		result[i] = derivative(x[i]);
	}
	return result;
}

double Yfunction::integral(const double x) {

	double y[YfunctionODEdim] = { 0 }; 

	YfunctionODEParams params;
	gsl_odeiv2_system sys = {Yfunction_rhs, NULL, YfunctionODEdim, &params};

	double from = -100.0;
	double h = std::min(1.e-6, tolerance);
	double to = x;
	while(from < to) {
		int status = gsl_odeiv2_evolve_apply (
			evolve, control, stepper,
			&sys,
			&from, to,
			&h, y
		);
		if (status != GSL_SUCCESS) {
			std::cerr << "ODE for Yfunction integration problems" << std::endl;
		}
	}
	double result = y[0];
	gsl_odeiv2_evolve_reset(evolve);
	gsl_odeiv2_step_reset(stepper);
    return result;
}