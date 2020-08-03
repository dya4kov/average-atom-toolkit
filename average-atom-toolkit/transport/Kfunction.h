#pragma once
#include <vector>
#include <cstddef>
#include <numeric-toolkit/interpolation/spline.h>

namespace aatk {
namespace transport {

using ::numtk::interpolation::Spline;

class Kfunction {
public:
	Kfunction(double n, double T, double M);
	~Kfunction();
	double operator()(double* e, double* tau, std::size_t npoints);
private:
	double  n, T, M;
	double  emin, emax;
	Spline* interpolation;
};

}
}
