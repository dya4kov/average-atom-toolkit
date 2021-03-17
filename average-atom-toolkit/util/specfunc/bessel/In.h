#pragma once
/**
* GNU Scientific Library In Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_In(const int n, const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct In {
	double operator() (const int n, const double x) {
		return gsl_sf_bessel_In(n, x);
	}
};

}
}
}
}