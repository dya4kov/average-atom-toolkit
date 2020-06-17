#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/semiclassic/atom.h>
#include <average-atom-toolkit/semiclassic/atom/ODE/action.h>

namespace aatk {
namespace semiclassic {

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using aatk::semiclassic::ODE::RHSAction;

double Atom::action(double energy, double lambda) {

	auto rpo = outerRP(energy, lambda);
	auto rpi = innerRP(energy, lambda);
	double xmax = std::sqrt(rpo[0]);
    double xmin = std::sqrt(rpi[0]);

    double result = 0.0;

    if (xmax <= xmin) return result;

    RHSAction rhs;
    Solver<PD853<RHSAction>> solver;
    solver.setTolerance(0.1*tolerance, 0.0);

    Array<RHSAction::dim> ay;
    ay[0] = rpo[1]; 
    ay[1] = rpo[2];
    ay[2] = 0.0;

    double l  = 0.5*lambda*lambda / r0 / r0;
    double e  = energy;

    rhs.set_e(e);
    rhs.set_l(l);
    rhs.set_eDens(densityInterpolation);
    
    solver.setStep(1e-10);
    solver.integrate(rhs, ay, xmax, xmin);

    result = 2.0*r0*std::sqrt(2.0)*ay[2];
    if (result < tolerance) result = 1e+10;

    return result;
}

}
}