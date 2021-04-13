#pragma once
/**
* GNU Scientific Library J0 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_J0(const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct J0 {
	double operator() (const double x) {
		return gsl_sf_bessel_J0(x);
	}
};

}
}
}
}