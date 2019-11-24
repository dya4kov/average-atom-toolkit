#pragma once
#include <numeric-toolkit/interpolation/base.h>

namespace numtk {
namespace interpolation {

class Polynomial : public Base {
public:
	Polynomial(
		const double*    x, 
		const double*    y, 
		std::size_t xysize,
		std::size_t  order
	);
	Polynomial(
		const std::vector<double>& x, 
		const std::vector<double>& y,
		std::size_t      order
	);
protected:
	double interpolate(std::size_t location, double x);
	double err;
};

}
}