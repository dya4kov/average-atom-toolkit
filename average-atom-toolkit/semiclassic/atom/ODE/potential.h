#pragma once
#include <cmath>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/interpolation/spline.h>

namespace aatk {
namespace semiclassic {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::interpolation::Spline;

struct RHSPotential {
	static const Dimension dim  = 2;

	RHSPotential() {}

    void set_eDens(Spline* _eDens) { eDens = _eDens; }

    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { 
        auto& density = *eDens;
        dydx[0] = 2.0*x*y[1];
        dydx[1] = x > 0 ? 2.0/x*density(x) : 0.0;
    }

private:
   
    Spline* eDens;

};

}
}
}