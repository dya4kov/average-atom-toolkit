#pragma once
#include <numeric-toolkit/interpolation/base.h>

namespace numtk {
namespace interpolation {

class Spline : public Base {
public:
	Spline(
		const double*    x, 
		const double*    y, 
		std::size_t xysize
	);
	Spline(
		const std::vector<double>& x, 
		const std::vector<double>& y
	);
protected:
	void   prepare();
	double interpolate(std::size_t location, double x);
	
	std::vector<double> y2;
};

}
}