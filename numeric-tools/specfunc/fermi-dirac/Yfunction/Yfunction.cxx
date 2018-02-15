#include <numeric-tools/specfunc/fermi-dirac/Yfunction.h>

using numtools::ODE::Array;
using numtools::ODE::Dimension;
using namespace numtools::specfunc;

Yfunction::Yfunction() : tolerance(1e-7) {
	solver.setTolerance(0.0, tolerance);
}

double Yfunction::operator()(const double& x) {
	return 0.5*FDhalf(x)*FDmhalf(x) + 1.5*integral(x);
}
double Yfunction::derivative(const double& x) {
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

void Yfunction::setTolerance(const double& eps) {
	tolerance = eps;
	solver.setTolerance(0.0, tolerance);
}

double Yfunction::integral(const double& x) {
	Array<ODE::RHSY::dim> y; y.fill(0.0);
    solver.integrate(rhs, y, -100.0, x);
    return y[ODE::RHSY::result];
}