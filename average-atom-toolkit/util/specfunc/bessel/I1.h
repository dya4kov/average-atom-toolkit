#pragma once
/**
* GNU Scientific Library I1 Bessel Function --- CXX Interface
*/
extern "C" { double gsl_sf_bessel_I1(const double x); }

namespace aatk {
namespace util {
namespace specfunc {
namespace bessel {

struct I1 {
	double operator() (const double x) {
		return gsl_sf_bessel_I1(x);
	}
};

}
}
}
}