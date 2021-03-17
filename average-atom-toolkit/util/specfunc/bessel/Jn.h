#pragma once
/**
* GNU Scientific Library Jn Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_Jn(const int n, const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct Jn {
	double operator() (const int n, const double x) {
		return gsl_sf_bessel_Jn(n, x);
	}
};

}
}
}
}