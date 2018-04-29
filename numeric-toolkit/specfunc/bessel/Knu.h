#pragma once
/**
* GNU Scientific Library Knu Bessel Function --- CXX Knuterface
*/
extern "C" { double gsl_sf_bessel_Knu(const double nu, const double x); }

namespace numtk {
namespace specfunc {
namespace bessel {

struct Knu {
	double operator() (const double& nu, const double& x) {
		return gsl_sf_bessel_Knu(nu, x);
	}
};

}
}
}