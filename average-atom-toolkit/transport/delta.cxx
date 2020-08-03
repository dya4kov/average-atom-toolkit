#include <cmath>
#include <cstddef>
#include <iostream>
#include <average-atom-toolkit/transport/delta.h>
#include <average-atom-toolkit/transport/bessel-spherical.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

namespace aatk {
namespace transport {

struct delta_EqParams {
    int l;
    double k;
    potential::Base* V;
    aatk::transport::Jl Jl;
	aatk::transport::Nl Nl;
};

int delta_odeFunction(double r, const double y[], double f[], void *params) {
    delta_EqParams* p = (delta_EqParams *)params;
    int    l = p->l;
    double k = p->k;
    auto&  V = *(p->V);
    auto& Jl = p->Jl;
    auto& Nl = p->Nl;
    double modulus = fabs(cos(y[0]) * Jl(l, k * r) - sin(y[0]) * Nl(l, k * r));
    f[0] = -1. / k * V(r) * modulus * modulus;
    return GSL_SUCCESS;
}

Delta::Delta(potential::Base& pot, double _rmax) : 
                    potential(pot),  rmax(_rmax) {}

double Delta::operator()(int l, double k) {
    const int dim = 1;
    double delta0 = potential.delta_eps(l, k);
    delta_EqParams params;
    params.l = l;
    params.k = k;
    params.V = &potential;

    const gsl_odeiv2_step_type * T
            = gsl_odeiv2_step_rk8pd;

    gsl_odeiv2_step * s
            = gsl_odeiv2_step_alloc (T, dim);
    gsl_odeiv2_control * c
            = gsl_odeiv2_control_y_new (1e-8, 0.0);
    gsl_odeiv2_evolve * e
            = gsl_odeiv2_evolve_alloc (dim);

    gsl_odeiv2_system sys = {delta_odeFunction, NULL, dim, &params};

    double from   = potential.r_eps();
	double to     = rmax;
    double h      = 1.e-6;
    double y[dim] = { delta0 };

    while(from < to) {
        int status = gsl_odeiv2_evolve_apply (e, c, s,
                                              &sys,
                                              &from, to,
                                              &h, y);

        if (status != GSL_SUCCESS) {
            std::cout << "ODE for delta integration problems" << std::endl;
        }
    }
    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
    return y[0];
}

}
}
