#pragma once
/**
* GNU Scientific Library K1 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_K1(const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct K1 {
	double operator() (const double x) {
		return gsl_sf_bessel_K1(x);
	}
};

}
}
}
}