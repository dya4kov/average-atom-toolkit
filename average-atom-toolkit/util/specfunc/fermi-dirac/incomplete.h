#pragma once
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

namespace aatk {
namespace util {
namespace specfunc {
/**
* Incomplete Fermi-Dirac functions
*/
template<typename Func>
struct FermiDiracInc {
	FermiDiracInc() : f() {}
	double operator() (const double x, const double y) {
		return f.value(x, y);
	}
private:
	Func f;
};

namespace FDI {

class HalfInc {
public:
	HalfInc(double tolerance = 1e-10);
	~HalfInc();
	double value(const double x, const double y);
private:	
	double tolerance;
	gsl_odeiv2_step    *stepper;
    gsl_odeiv2_control *control;
    gsl_odeiv2_evolve   *evolve;
};

class ThreeHalfInc {
public:
	ThreeHalfInc(double tolerance = 1e-10);
	~ThreeHalfInc();
	double value(const double x, const double y);
private:	
	double tolerance;
	gsl_odeiv2_step    *stepper;
    gsl_odeiv2_control *control;
    gsl_odeiv2_evolve   *evolve;
};

}
}
}
}