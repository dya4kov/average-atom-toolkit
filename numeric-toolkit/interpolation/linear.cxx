#include <numeric-toolkit/interpolation/linear.h>

using namespace numtk::interpolation;

Linear::Linear(
	const double*         x, 
	const double*         y, 
	std::size_t      xysize
) : Base(x, y, xysize, 2) { }

Linear::Linear(
	const std::vector<double>&  x, 
	const std::vector<double>&  y
) : Base(x, y, 2) { }

double Linear::interpolate(std::size_t ix, double x) {
	if (xx[ix] == xx[ix + 1]) return yy[ix];
	else return yy[ix] + ((x - xx[ix])/(xx[ix + 1] - xx[ix]))*(yy[ix + 1] - yy[ix]);
}
