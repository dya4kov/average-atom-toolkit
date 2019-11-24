#include <numeric-toolkit/interpolation/spline.h>
#include <cmath>

using namespace numtk::interpolation;

Spline::Spline(
	const double*    x, 
	const double*    y, 
	std::size_t xysize
) : Base(x, y, xysize, 2), y2(xysize) { prepare(); }

Spline::Spline(
	const std::vector<double>&  x, 
	const std::vector<double>&  y
) : Base(x, y, 2), y2(x.size()) { prepare(); }

double Spline::interpolate(std::size_t ix, double x) {
	std::size_t ilow = ix; 
	std::size_t iup  = ix + 1;
	double y, h, b, a;
	h = xx[iup] - xx[ilow];
	if (h == 0.0) throw("Bad input to spline interpolation");
	a = (xx[iup] - x)/h;
	b = (x - xx[ilow])/h;
	y = a*yy[ilow] + b*yy[iup] 
	  + ((a*a*a - a)*y2[ilow]
	  +	 (b*b*b - b)*y2[iup])
	  *  (h*h)/6.0;
	return y;
}

void Spline::prepare() {
	double p, qn, sig, un;
	std::size_t n = y2.size();
	std::vector<double> u(n - 1);
	y2[0] = u[0] = 0.0;
	for (std::size_t i = 1; i < n - 1; ++i) {
		sig   = (xx[i] - xx[i - 1])/(xx[i + 1] - xx[i - 1]);
		p     = sig*y2[i - 1] + 2.0;
		y2[i] = (sig - 1.0)/p;
		u[i]  = (yy[i + 1] - yy[i])/(xx[i + 1] - xx[i]) - (yy[i] - yy[i - 1])/(xx[i] - xx[i - 1]);
		u[i]  = (6.0*u[i]/(xx[i + 1] - xx[i - 1]) - sig*u[i - 1])/p;
	}
	y2[n - 1] = 0.0;
	for (std::size_t k = n - 2; k > 0; --k) {
		y2[k] = y2[k]*y2[k + 1] + u[k];
	}
	y2[0] = y2[0]*y2[1] + u[0];
}
