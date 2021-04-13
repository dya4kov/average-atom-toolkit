#pragma once
/**
* GNU Scientific Library Yn Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_Yn(const int n, const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct Yn {
	double operator() (const int n, const double x) {
		return gsl_sf_bessel_Yn(n, x);
	}
};

}
}
}
}