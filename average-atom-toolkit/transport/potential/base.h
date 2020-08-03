#pragma once
#include <vector>
#include <cstddef>

#include <numeric-toolkit/specfunc/gamma/double-fact.h>
#include <average-atom-toolkit/transport/bessel-spherical.h>

namespace aatk {
namespace transport {
namespace potential {

class Base {
public:
	Base(double _eps = 1.e-6) : eps(_eps) {}
	double r_eps() { return eps; }

	virtual double  delta_eps(int l, double k) { return 0.0; }
	virtual double  operator()(double r) { return 0.0; }
	        double* operator()(double* r, double* V, std::size_t npoints);
	virtual ~Base() {}
protected:

	// initial coordinate for differential equation (5)
	double eps;
	// special functions
	::numtk::specfunc::gamma::DoubleFact doubleFact;
	aatk::transport::Jl Jl;
	aatk::transport::Nl Nl;
};

}
}
}