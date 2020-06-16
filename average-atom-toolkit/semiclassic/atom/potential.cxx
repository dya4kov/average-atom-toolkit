#include <cmath>
#include <numeric>
#include <algorithm>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/semiclassic/atom.h>
#include <average-atom-toolkit/semiclassic/atom/ODE/potential.h>

namespace aatk {
namespace semiclassic {

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using aatk::semiclassic::ODE::RHSPotential;

std::vector<double> Atom::potential(const std::vector<double>& x) {
	return std::vector<double>(x.size(), 0.0);
}
double Atom::potential(double x) {
	return 0.0;
}
void Atom::potential(const double* x, double* result, std::size_t n) {
	std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
       [x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    RHSPotential rhs;

    rhs.set_eDens(densityInterpolation);

    Array<RHSPotential::dim> phi; 
    double xFrom = 1.0;
    phi.fill(0.0);

    Solver<PD853<RHSPotential>> solver;
    solver.setTolerance(0.0, 0.1*tolerance);

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        solver.setStep(tolerance);
        solver.integrate(rhs, phi, xFrom, xTo);
        result[i] = phi[0]/x[i];
        xFrom = xTo;
    }
}

}
}