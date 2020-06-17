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

struct RHSAction {
	static const Dimension dim    = 3;
    static const Dimension result = 2;

	RHSAction() : energyArg(1.0), lambdaArg(1.0) {}

    void set_l (const double& _l)  { lambdaArg = _l; }
    void set_e (const double& _e)  { energyArg = _e; }
    void set_eDens(Spline* _eDens) { eDens = _eDens; }
		
	void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { 
        auto& density = *eDens;
        dydx[0] = 2.0*x*y[1];
        dydx[1] = x > 0 ? 2.0/x*density(x) : 0.0;
        double p = energyArg*x*x + y[0] - lambdaArg/(x*x);
        dydx[2] = p > 0.0 ? -std::sqrt(p) : 0.0;
    }
private:

	double energyArg;
	double lambdaArg;
    Spline* eDens;
};

}
}
}