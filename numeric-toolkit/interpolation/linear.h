#pragma once
#include <numeric-toolkit/interpolation/base.h>

namespace numtk {
namespace interpolation {

class Linear : public Base {
public:
	Linear(
		const double*    x, 
		const double*    y, 
		std::size_t xysize
	);
	Linear(
		const std::vector<double>& x, 
		const std::vector<double>& y
	);
protected:
	double interpolate(std::size_t location, double x); 
};

}
}