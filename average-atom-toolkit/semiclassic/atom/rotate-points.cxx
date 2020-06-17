#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/semiclassic/atom.h>
#include <average-atom-toolkit/semiclassic/atom/ODE/rotate-points.h>

namespace aatk {
namespace semiclassic {

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;

using aatk::semiclassic::ODE::RHSRPouter;
using aatk::semiclassic::ODE::RHSRPinner;

std::array<double, 3> Atom::innerRP(double energy, double lambda) {

	std::array<double, 3> result = {1.0, 0, 0};

    double l = 0.5*lambda*lambda / r0 / r0;
    double e = energy;

    RHSRPinner rhsRP;
    rhsRP.reset();
    rhsRP.set_e(e);
    rhsRP.set_l(l);
    rhsRP.set_eDens(densityInterpolation);

    double t = 1e-10;
    double h = 1e-9;

    auto rpo = outerRP(energy, lambda); 

    double xmin = 0.0;
    double xmax = std::sqrt(rpo[0]);

    Array<RHSRPinner::dim> RPy;
    RPy[0] = rpo[1];
    RPy[1] = rpo[2];

    double error = std::abs(xmax - xmin);
    Solver<PD853<RHSRPinner>> solverRP;
    int nStep = 0;

    while (error > 0.1*tolerance && nStep < 20) {

        solverRP.setStep(h);
        solverRP.setTolerance(0.0, t);
        solverRP.integrate(rhsRP, RPy, xmax, xmin);

        xmin = rhsRP.xDown();
        xmax = rhsRP.xUp();

        error = std::abs(xmax - xmin);
        RPy   = rhsRP.yUp();
        h     = error / 11.0;
        t     = h / 21.0;

        ++nStep;
    }

    result[0] = xmax*xmax;
    result[1] = rhsRP.yUp()[0];
    result[2] = rhsRP.yUp()[1];

	return result;
}

std::array<double, 3> Atom::outerRP(double energy, double lambda) {

	std::array<double, 3> result = {1.0, 0, 0};

    double l = 0.5*lambda*lambda / r0 / r0;
    double e = energy;

    if (e > l) return result;

    RHSRPouter rhsRP;
    rhsRP.reset();
    rhsRP.set_e(e);
    rhsRP.set_l(l);
    rhsRP.set_eDens(densityInterpolation);

    double t = 1e-10;
    double h = 1e-9;

    double xmin = 0.0;
    double xmax = 1.0;

    Array<RHSRPouter::dim> RPy; RPy.fill(0.0);

    double error = std::abs(xmax - xmin);
    Solver<PD853<RHSRPouter>> solverRP;
    int nStep = 0;

    while (error > 0.1*tolerance && nStep < 20) {

        solverRP.setStep(h);
        solverRP.setTolerance(0.0, t);
        solverRP.integrate(rhsRP, RPy, xmax, xmin);

        xmin = rhsRP.xDown();
        xmax = rhsRP.xUp();

        error = std::abs(xmax - xmin);
        RPy   = rhsRP.yUp();
        h     = error  / 11.0;
        t     = h / 21.0;

        ++nStep;
    }

    result[0] = xmin*xmin;
    result[1] = rhsRP.yDown()[0];
    result[2] = rhsRP.yDown()[1];
	return result;
}

}
}