#pragma once
/**
* GNU Scientific Library K0 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_K0(const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct K0 {
	double operator() (const double x) {
		return gsl_sf_bessel_K0(x);
	}
};

}
}
}
}