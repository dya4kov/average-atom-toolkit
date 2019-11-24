#pragma once
#include <vector>
#include <cstddef>

namespace numtk {
namespace interpolation {

class Base {
public:
	Base(
		const double*        x, 
		const double*        y, 
		std::size_t     xysize, 
		std::size_t interpSize
	);
	Base(
		const std::vector<double>& x, 
		const std::vector<double>& y, 
		std::size_t interpSize
	);
	Base(const Base& interpolator);
	Base& operator=(const Base& interpolator);

	double operator()(const double x);

protected:

	std::size_t locate(const double x);
	std::size_t hunt(const double x);

	virtual double interpolate(std::size_t location, double x); 

	bool        correlated;
	std::size_t corDistance;
	std::size_t interpSize;
	std::size_t isaved;

	std::vector<double> xx, yy;

};

}
}
