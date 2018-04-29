#pragma once
/**
* GNU Scientific Library Y1 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_Y1(const double x); }

namespace numtk {
namespace specfunc {
namespace bessel {

struct Y1 {
	double operator() (const double& x) {
		return gsl_sf_bessel_Y1(x);
	}
};

}
}
}