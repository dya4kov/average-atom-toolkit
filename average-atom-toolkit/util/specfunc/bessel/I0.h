#pragma once
/**
* GNU Scientific Library I0 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_I0(const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct I0 {
	double operator() (const double x) {
		return gsl_sf_bessel_I0(x);
	}
};

}
}
}
}