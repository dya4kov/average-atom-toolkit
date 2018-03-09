#pragma once

#include <numeric-toolkit/ODE/types.h>

namespace numtk {
namespace ODE {
namespace stepper {

// Dormand-Prince fifth-order Runge-Kutta step with monitoring of 
// local truncation error to ensure accuracy and adjust stepsize.

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

struct PD5_constants {
	const double
	c2,   c3,  c4,  c5, 
	a21, a31, a32, a41, a42, a43, 
	a51, a52, a53, a54, 
	a61, a62, a63, a64, a65, 
	a71, a73, a74, a75, a76, 
	e1, e3, e4, e5, e6, e7;
	PD5_constants() :
	 c2( 0.2),
	 c3( 0.3),
	 c4( 0.8),
	 c5( 8.0/9.0),
	a21( 0.2),
	a31( 3.0/40.0),
	a32( 9.0/40.0),
	a41( 44.0/45.0),
	a42(-56.0/15.0),
	a43( 32.0/9.0),
	a51( 19372.0/6561.0),
	a52(-25360.0/2187.0),
	a53( 64448.0/6561.0),
	a54(-212.0/729.0),
	a61( 9017.0/3168.0),
	a62(-355.0/33.0),
	a63( 46732.0/5247.0),
	a64( 49.0/176.0),
	a65(-5103.0/18656.0),
	a71( 35.0/384.0),
	a73( 500.0/1113.0),
	a74( 125.0/192.0),
	a75(-2187.0/6784.0),
	a76( 11.0/84.0),
	 e1( 71.0/57600.0),
	 e3(-71.0/16695.0),
	 e4( 71.0/1920.0),
	 e5(-17253.0/339200.0),
	 e6( 22.0/525.0),
	 e7(-1.0/40.0) {}
};


template<typename RHStype>
class PD5 : public PD5_constants {
public:
	typedef RHStype RHS;

	PD5(
		Array<RHStype::dim>  &_y, 
		Array<RHStype::dim>  &_dydx, 
		double               &_x, 
		double               &_tolAbs, 
		double               &_tolRel, 
		double               &_EPS) : 
    	y(_y), dydx(_dydx), x(_x), 
    	tolAbs(_tolAbs), tolRel(_tolRel), EPS(_EPS),
    	PD5_constants() {}

	void dy(const double h, RHStype &rhs);
	
	void   makeStep(const double hTry, RHStype &rhs);
	double lastStep() { return hDid;  }
	double nextStep() { return hNext; }

	double error();

	struct Controller {
		Controller();
		bool success(const double err, double &h);
		double hNext, errOld;
		bool reject;
	};
private:
	double               &x;
	Array<RHStype::dim>  &y, &dydx;
	Array<RHStype::dim>  k2,k3,k4,k5,k6;
	Array<RHStype::dim>  dydxnew;
	Array<RHStype::dim>  yOut, yErr;
	Controller           controller;

	double &tolAbs, &tolRel, &EPS;
	double hDid, hNext;
};

template<class RHStype>
void PD5<RHStype>::makeStep(const double hTry, RHStype &rhs) {
	double h = hTry; // Set stepsize to the initial trial value
	for (;;) {
		dy(h, rhs); // Take a step
		double err = error(); // Evaluate accuracy
		if (controller.success(err, h)) break;
		// Else step rejected. Try again with
        // reduced h set by controller
		if (std::abs(h) <= EPS) break;
	}
	dydx = dydxnew; // Reuse last derivative evaluation for next step
	y = yOut;
	hDid = h;
	x += hDid;
	hNext = controller.hNext;
}

template<class RHStype>
void PD5<RHStype>::dy(const double h, RHStype &rhs) {
	Array<RHStype::dim> yTemp; yTemp.fill(0.0);
	Dimension i;
	for (i = 0; i < RHStype::dim; ++i) { // First step
		yTemp[i] = y[i] + h*a21*dydx[i];
	}
	rhs(x + c2*h, yTemp, k2); // Second step
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a31*dydx[i] + a32*k2[i]);
	}
	rhs(x + c3*h, yTemp, k3); // Third step
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a41*dydx[i] + a42*k2[i] + a43*k3[i]);
	}
	rhs(x + c4*h, yTemp, k4); // Fourth step
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a51*dydx[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
	}
	rhs(x + c5*h, yTemp, k5); // Fifth step
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a61*dydx[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
	}
	double xph = x + h;
	rhs(xph, yTemp, k6); // Sixth step
	for (i = 0; i < RHStype::dim; ++i) {
		// Accumulate increments with proper weights
		yOut[i] = y[i] + h*(a71*dydx[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
	}
	rhs(xph,yOut,dydxnew); // Will also be first evaluation for next step.
	//Estimate error as difference between fourth- and fifth-order methods.
	for (i = 0; i < RHStype::dim; ++i) {
		yErr[i]=h*(e1*dydx[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*dydxnew[i]);
	}
}

template<class RHStype>
double PD5<RHStype>::error() {
	// Use yerr to compute norm of scaled error estimate.
	// A value less than one means the step wassuccessful
	double err = 0.0;
	double sk;
	for (Dimension i = 0; i < RHStype::dim; ++i) {
		sk = tolAbs + tolRel*std::max(std::abs(y[i]),std::abs(yOut[i]));
		sk = yErr[i]/sk;
		sk *= sk;
		err += sk;
	}
	return std::sqrt(err/RHStype::dim);
}

template<class RHStype>
PD5<RHStype>::Controller::Controller() : reject(false), errOld(1.0e-4) {}

template <class RHStype>
bool PD5<RHStype>::Controller::success(const double err, double &h) {
	// Returns true if err <= 1, false otherwise. If step was successful, 
	// sets hNext to the estimated optimal stepsize for the next step.
	// If the step failed, reduces h appropriately for another try.
	const double 
	beta = 0.04, 
	alpha = 0.2 - beta*0.75, 
	safe = 0.9, 
	minscale = 0.2,
	maxscale = 10.0;
	// Set beta to a nonzero value for PI control. beta = 0.04-0.08 is a good default.
	double scale;
	if (err <= 1.0) { // Step succeeded. Compute hNext.
		if (err == 0.0) {
			scale = maxscale;
		}
		else { // PI control if beta != 0.
			scale = safe*std::pow(err, -alpha)*std::pow(errOld, beta);
			// Ensure minscale <= hNext/h <= maxscale.
			if (scale < minscale) scale = minscale;  
			if (scale > maxscale) scale = maxscale;
		}
		// Don’t let step increase if last one was rejected
		if (reject) {
			hNext = h*std::min(scale,1.0);
		}
		else {
			hNext = h*scale;
		}
		errOld = std::max(err, 1.0e-4); // Bookkeeping for next call
		reject = false;
		return true;
	}
	else {
		scale = std::max(safe*std::pow(err, -alpha), minscale);
		h *= scale;
		reject = true;
		return false;
	}
}

} // end of namespace stepper
} // end of namespace ODE
} // end of namespace numtk