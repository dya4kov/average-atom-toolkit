#pragma once
/**
* GNU Scientific Library double factorial n!! --- CXX Interface
*/
extern "C" { double gsl_sf_doublefact(const unsigned int n); }

namespace aatk {
namespace util {
namespace specfunc {
namespace gamma {

struct DoubleFact {
	double operator() (const unsigned int& n) {
		return gsl_sf_doublefact(n);
	}
};

}
}
}
}