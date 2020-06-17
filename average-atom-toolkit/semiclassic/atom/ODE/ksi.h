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

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

struct RHSksi {
	static const Dimension dim    = 3;
    static const Dimension result = 2;

	RHSksi() : energyArg(1.0), lambdaArg(1.0) {}

    void set_eDens(Spline* _eDens) { eDens = _eDens; }
    void set_l (const double& _l)  { lambdaArg = _l; }
    void set_e (const double& _e)  { energyArg = _e; }
		
	void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { 
        auto&  density = *eDens;
        double x2 = x*x;
        double p = std::abs(energyArg*x2 + y[0] - lambdaArg/x2);

        dydx[0] = 2.0*x*y[1];
        dydx[1] = x > 0 ? 2.0*density(x)/x : 0.0;
        dydx[2] = -std::sqrt(p);
    }
private:

	double energyArg;
	double lambdaArg;
    Spline* eDens;
};

}
}
}