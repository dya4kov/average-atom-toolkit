#include <cmath>
#include <average-atom-toolkit/atom/base.h>

namespace aatk {
namespace atom {

// Atom::Atom(ConstructorArgs args) {
// 	reset(args);
// }

// old style
Atom::Atom(
	double _V,
	double _T,
	double _Z,
	double _tolerance,
	int    _meshSize
) {
	reset(_V, _T, _Z, _tolerance, _meshSize);
}

// void Atom::reset(ConstructorArgs args) {
// 	reset(
// 		args.V, 
// 		args.T, 
// 		args.Z, 
// 		args.tolerance, 
// 		args.meshSize
// 	);
// }

void Atom::reset(
	double _V,
	double _T,
	double _Z,
	double _tolerance,
	int    _meshSize
) {
	V = _V;
	T = _T;
	Z = _Z;
	M = 0.0;
	tolerance = _tolerance;
	meshSize = _meshSize;
	r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
}

Atom::~Atom() {}

double Atom::radius() { return r0; } 

double Atom::volume() { return V; }

double Atom::temperature() { return T;}

double Atom::Znucleus() { return Z; }

double Atom::ZfreeElectrons() { return 0.0; }

double Atom::chemicalPotential() { return M; }

void Atom::U(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i] = 0.0;
	}
	return;
}

std::vector<double> Atom::U(const std::vector<double>& x) {
	std::vector<double> y(x.size());
	U(x.data(), y.data(), x.size());
	return y;
}

double Atom::U(double x) {
	double y;
	U(&x, &y, 1);
	return y;
}

void Atom::xU(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i] = 0.0;
	}
	return;
}

std::vector<double> Atom::xU(const std::vector<double>& x) {
	std::vector<double> y(x.size());
	xU(x.data(), y.data(), x.size());
	return y;
}

double Atom::xU(double x) {
	double y;
	xU(&x, &y, 1);
	return y;
}

void Atom::x2dU(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i]  = 0.0;
	}
	return;
}

std::vector<double> Atom::x2dU(const std::vector<double>& x) {
	std::vector<double> y(x.size());
	x2dU(x.data(), y.data(), x.size());
	return y;
}

double Atom::x2dU(double x) {
	double y;
	x2dU(&x, &y, 1);
	return y;
}

void Atom::electronDensity(const double* x, double* dens, std::size_t n, double eb) {
	for (std::size_t i = 0; i < n; ++i) {
		dens[i] = 0.0;
	}
	return;
}

std::vector<double> Atom::electronDensity(const std::vector<double>& x, double eb) {
	std::vector<double> dens(x.size());
	electronDensity(x.data(), dens.data(), x.size(), eb);
	return dens;
}

double Atom::electronDensity(double x, double eb) {
	double dens;
	electronDensity(&x, &dens, 1, eb);
	return dens;
}

std::vector<double> Atom::waveFunction(
	const std::vector<double>& x, 
	double                     e, 
	double                     lambda
) {
	std::vector<double> dens(x.size());
	waveFunction(x.data(), dens.data(), x.size(), e, lambda);
	return dens;
}

double Atom::waveFunction(
	double x, 
	double e, 
	double lambda
) {
	double wf;
	waveFunction(&x, &wf, 1, e, lambda);
	return wf;
}

void Atom::waveFunction(
	const double *x, 
	double       *y, 
	std::size_t   n, 
	double        e, 
	double        lambda
) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i] = 0.0;
	}
	return;
}

double Atom::innerRP(double e, double lambda) {
	return 0.0;
}
double Atom::outerRP(double e, double lambda) {
	return 0.0;
}
double Atom::action(double e, double lambda) {
	return 0.0;
}
double Atom::energyLevel(int n, int l) {
	return 0.0;
}
double Atom::electronStatesDiscrete(int n, int l) {
	return 0.0;
}
double Atom::electronStatesDiscrete(int n) {
	return 0.0;
}
double Atom::electronStatesDiscrete() {
	return 0.0;
}
double Atom::electronStatesContinuous() {
	return 0.0;
}

}
}