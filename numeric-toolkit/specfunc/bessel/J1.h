#pragma once
/**
* GNU Scientific Library J1 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_J1(const double x); }

namespace numtk {
namespace specfunc {
namespace bessel {

struct J1 {
	double operator() (const double& x) {
		return gsl_sf_bessel_J1(x);
	}
};

}
}
}