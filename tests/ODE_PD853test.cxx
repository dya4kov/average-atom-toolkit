#include <iostream>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

using numtk::ODE::Array;
using numtk::ODE::Dimension;

struct rhs_harmonic {
	double omega;
	static const Dimension dim = 2;
	rhs_harmonic(double _omega) : omega(_omega) {}
	void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) {
		dydx[0] = y[1];
		dydx[1] = -omega*omega*y[0];
	}
};

int main() {
	numtk::ODE::Solver<numtk::ODE::stepper::PD853<rhs_harmonic>> ode;
	rhs_harmonic rhs(1.0);
	Array<rhs_harmonic::dim> y;
	y[0] = 1.0;
	y[1] = 0.0;
	//auto y = ode.Integrate(rhs, y, 0.0, 3.14);
	for (double x = 0.1; x < 6.35; x += 0.1) {
		ode.integrate(rhs, y, x - 0.1, x);
		std::cout <<  x << "     " << y[0] << "    " << std::cos(x) << std::endl;
	}
}