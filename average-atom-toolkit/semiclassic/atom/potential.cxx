#include <cmath>
#include <numeric>
#include <algorithm>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

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

struct potentialODEParams {
    Spline* eDens;
};

int potentialODE(double x, const double y[], double dydx[], void *params) {
    auto  p = (potentialODEParams *) params;
    auto& density = *(p->eDens);

    dydx[0] = 2.0*x*y[1];
    dydx[1] = x > 0 ? 2.0/x*density(x) : 0.0;

    return GSL_SUCCESS;
}

std::vector<double> Atom::potential(const std::vector<double>& x) {
    std::vector<double> result(x.size(), 0.0);
    potential(x.data(), result.data(), x.size());
	return result;
}

double Atom::potential(double x) {
    double result;
    potential(&x, &result, 1);
	return result;
}
void Atom::potential(const double* x, double* result, std::size_t n) {
	std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
       [x](std::size_t i1, std::size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    // RHSPotential rhs;

    // rhs.set_eDens(densityInterpolation);

    // Array<RHSPotential::dim> phi; 
    // double xFrom = 1.0;
    // phi.fill(0.0);

    // Solver<PD853<RHSPotential>> solver;
    // solver.setTolerance(0.0, 0.1*tolerance);

    // for (auto i : idx) {
    //     double xTo = std::sqrt(x[i]);
    //     solver.setStep(tolerance);
    //     solver.integrate(rhs, phi, xFrom, xTo);
    //     result[i] = phi[0]/x[i];
    //     xFrom = xTo;
    // }

    const int dim = 2;
    potentialODEParams params;
    params.eDens = densityInterpolation;

    const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rk8pd;

    gsl_odeiv2_step*    stepper = gsl_odeiv2_step_alloc (stepType, dim);
    gsl_odeiv2_control* control = gsl_odeiv2_control_y_new (0.0, 0.1*tolerance);
    gsl_odeiv2_evolve*   evolve = gsl_odeiv2_evolve_alloc (dim);

    gsl_odeiv2_system sys = {potentialODE, NULL, dim, &params};

    double from     = 1.0;
    double h        = -1.e-6;
    double phi[dim] = { 0 };

    for (auto i : idx) {
        double to = std::sqrt(x[i]);
        while(from > to) {
            int status = gsl_odeiv2_evolve_apply (evolve, control, stepper,
                                                  &sys,
                                                  &from, to,
                                                  &h, phi);
            if (status != GSL_SUCCESS) {
                std::cout << "ODE for potential integration problems" << std::endl;
            }
        }
        result[i] = phi[0]/x[i];
        from = to;
    }

    gsl_odeiv2_evolve_free (evolve);
    gsl_odeiv2_control_free (control);
    gsl_odeiv2_step_free (stepper);
}

}
}