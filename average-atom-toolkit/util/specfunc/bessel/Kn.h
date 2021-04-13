#pragma once
/**
* GNU Scientific Library Kn Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_Kn(const int n, const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct Kn {
	double operator() (const int n, const double x) {
		return gsl_sf_bessel_Kn(n, x);
	}
};

}
}
}
}