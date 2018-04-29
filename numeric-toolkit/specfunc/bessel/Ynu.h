#pragma once
/**
* GNU Scientific Library Ynu Bessel Function --- CXX Ynuterface
*/
extern "C" { double gsl_sf_bessel_Ynu(const double nu, const double x); }

namespace numtk {
namespace specfunc {
namespace bessel {

struct Ynu {
	double operator() (const double& nu, const double& x) {
		return gsl_sf_bessel_Ynu(nu, x);
	}
};

}
}
}