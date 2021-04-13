#pragma once
/**
* GNU Scientific Library Y0 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_Y0(const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct Y0 {
	double operator() (const double x) {
		return gsl_sf_bessel_Y0(x);
	}
};

}
}
}
}