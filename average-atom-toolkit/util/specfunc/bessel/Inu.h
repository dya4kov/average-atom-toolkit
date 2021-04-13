#pragma once
/**
* GNU Scientific Library Inu Bessel Function --- CXX Inuterface
*/
extern "C" { double gsl_sf_bessel_Inu(const double nu, const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct Inu {
	double operator() (const double nu, const double x) {
		return gsl_sf_bessel_Inu(nu, x);
	}
};

}
}
}
}