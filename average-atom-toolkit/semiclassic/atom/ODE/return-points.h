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

struct RHSRPouter {
	static const Dimension dim = 2;

	double     xUp()   { return _xUp;   }
	double     xDown() { return _xDown; }
	Array<dim> yUp()   { return _yUp;   }
	Array<dim> yDown() { return _yDown; }

	RHSRPouter() : energyArg(1.0), lambdaArg(1.0) {
        _yOld.fill(-1.0); p1 = -1.0;
        _xOld =     1.0;  p2 = -1.0;
    }


    void reset  () { _xUp = 0.0; _xDown = 0.0; }
    void set_l  (const double& _l) { lambdaArg = _l; }
    void set_e  (const double& _e) { energyArg = _e; }
    void set_eDens(Spline* _eDens) { eDens = _eDens; }

	void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { 
        auto& density = *eDens;

        dydx[0] = 2.0*x*y[1];
        dydx[1] = x > 0 ? 2.0/x*density(x) : 0.0;
        
        p2 = energyArg*x*x + y[0] - lambdaArg / (x*x);
        
        if (p1 < 0 && p2 >= 0 && x < _xOld) {
            _yUp = _yOld; _yDown = y;
            _xUp = _xOld; _xDown = x;
        }

        p1 = p2;
        _yOld = y;
        _xOld = x;
    }
private:

	double _xOld;
    double _xUp;
	double _xDown;

	Array<dim> _yOld;
	Array<dim> _yUp;
	Array<dim> _yDown;

	double energyArg;
	double lambdaArg;
	double p1;
	double p2;

    Spline* eDens;
};

struct RHSRPinner {

	static const Dimension dim = 2;

    double     xUp()   { return _xUp;   }
	double     xDown() { return _xDown; }
	Array<dim> yUp()   { return _yUp;   }
	Array<dim> yDown() { return _yDown; }
	
	RHSRPinner() : 
        energyArg(1.0), 
        lambdaArg(1.0) 
    {
        _yOld.fill(-1.0); 
        p1   = -1.0;
		p2   = -1.0;
    }

    void reset  () { _xUp = 0.0; _xDown = 0.0; }
    void set_l  (const double& _l) { lambdaArg = _l; }
    void set_e  (const double& _e) { energyArg = _e; }
    void set_eDens(Spline* _eDens) { eDens = _eDens; }

	void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) {
        auto& density = *eDens;
        
        dydx[0] = 2.0*x*y[1];
        dydx[1] = x > 0 ? 2.0/x*density(x) : 0.0;

        p2 = energyArg*x*x + y[0] - lambdaArg / (x*x);

        if (p1 >= 0 && p2 < 0 && x < _xOld) {
            _yUp = _yOld; _yDown = y;
            _xUp = _xOld; _xDown = x;
        }

        p1 = p2;
        _yOld = y;
        _xOld = x; 
    }
private:

	double _xOld;
    double _xUp;
	double _xDown;

	Array<dim> _yOld;
	Array<dim> _yUp;
	Array<dim> _yDown;

	double energyArg;
	double lambdaArg;
	double p1;
	double p2;

    Spline* eDens;
};

}
}
}