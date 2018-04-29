#pragma once
#include <cmath>
#include <limits>
#include <iostream>
#include <array>

#include <numeric-toolkit/ODE/types.h>

namespace numtk {
namespace ODE {

template<typename Stepper, Dimension dim = Stepper::RHS::dim>
class Solver {
public:
	Solver(double _hStart = 1e-6, double _hMin = std::numeric_limits<double>::epsilon(), unsigned _maxStep = 50000) : 
	tolAbs(1e-6), tolRel(0.0), EPS(0.0),
	stepper(y, dydx, x, tolAbs, tolRel, EPS) {
		y.fill(0.0); dydx.fill(0.0); x = 0.0;
		hStart  = _hStart;
		hMin    = _hMin;
		maxStep = _maxStep;
		nOk = 0;
		nBad = 0;
		EPS = std::numeric_limits<double>::epsilon();
	}
	void setTolerance(const double _tolAbs, const double _tolRel) {
		tolAbs = _tolAbs;
		tolRel = _tolRel;
	}	
	void setStep(const double _hStart) { hStart = _hStart; }
	void integrate(typename Stepper::RHS &rhs, Array<dim> &_yStart, 
		 		   const double& _xFrom, const double& _xTo) 
	{
		yStart = _yStart;
		y      =  yStart;
		xFrom  = _xFrom;
		xTo    = _xTo;
		x      =  xFrom;

		h = std::copysign(std::abs(hStart), xTo - xFrom);

		rhs(x, y, dydx); // store initial values

		for (nStep = 0; nStep < maxStep; ++nStep) {
			if ((x + h*1.0001 - xTo)*(xTo - xFrom) > 0.0) {
				h = xTo - x; // if stepsize can overshoot, decrease
			}
			stepper.makeStep(h, rhs);
			if (stepper.lastStep() == h) ++nOk; else ++nBad;
			if ((x - xTo)*(xTo - xFrom) >= 0.0) { // Are we done?
				_yStart = y; // update ystart
				return; // normal exit
			}
			//if (std::abs(stepper.nextStep()) <= hMin) 
			//	std::cerr << "Step size too small in ODEsolver" << std::endl;
			h = stepper.nextStep();
		}
		// std::cerr << "Too many steps in routine ODEsolver" << std::endl;
	}

private:
	unsigned maxStep;
	unsigned nOk;
	unsigned nBad;
	unsigned nStep;

	double tolAbs, tolRel;
	double xFrom,  xTo, x;
	double hStart, hMin, h;
	double EPS;

	Array<dim> y, dydx, yStart;
	Stepper stepper;

};

}
}