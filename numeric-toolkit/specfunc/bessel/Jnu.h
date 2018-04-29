#pragma once
/**
* GNU Scientific Library Jnu Bessel Function --- CXX Jnuterface
*/
extern "C" { double gsl_sf_bessel_Jnu(const double nu, const double x); }

namespace numtk {
namespace specfunc {
namespace bessel {

struct Jnu {
	double operator() (const double& nu, const double& x) {
		return gsl_sf_bessel_Jnu(nu, x);
	}
};

}
}
}