#include <numeric-toolkit/interpolation/polynomial.h>
#include <cmath>

using namespace numtk::interpolation;

Polynomial::Polynomial(
	const double*    x, 
	const double*    y, 
	std::size_t xysize,
	std::size_t  order
) : Base(x, y, xysize, order + 1) { }

Polynomial::Polynomial(
	const std::vector<double>&  x, 
	const std::vector<double>&  y,
	std::size_t       order
) : Base(x, y, order + 1) { }

double Polynomial::interpolate(std::size_t ix, double x) {
	double y, den, dif, dift, ho, hp, w;
	const double* xa = xx.data() + ix;
	const double* ya = yy.data() + ix;
	std::vector<double> c(interpSize), d(interpSize);
	dif = std::abs(x - xa[0]);
	// ns - the closest to x table entry
	std::size_t ns = 0;
	for (std::size_t i = 0; i < interpSize; ++i) {
		dift = std::abs(x - xa[i]);
		if (dift < dif) { ns = i; dif = dift; }
		c[i] = ya[i];
		d[i] = ya[i];
	}
	y = ya[ns--]; // initial approximation to y
	for (std::size_t m = 1; m < interpSize; ++m) {
		for (std::size_t i = 0; i < interpSize - m; ++i) {
			ho = xa[i] - x;
			hp = xa[i + m] - x;
			w = c[i + 1] - d[i];
			den = ho - hp;
			if (den == 0.0) throw("Polynomial interpolation error");
			den = w/den;
			d[i] = hp*den;
			c[i] = ho*den;
		}
		err = 2*(ns + 1) < interpSize - m ? c[ns + 1] : d[ns--];
		y += err;
	}
	return y;
}
