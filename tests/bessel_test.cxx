#include <numeric-toolkit/specfunc/bessel/Jnu.h>
#include <numeric-toolkit/specfunc/bessel/Ynu.h>
#include <numeric-toolkit/specfunc/bessel/Knu.h>
#include <iostream>
#include <cmath>


int main() {
	double gamma = 2.6789385347077476336557;
	numtk::specfunc::bessel::Jnu Jnu;
	numtk::specfunc::bessel::Ynu Ynu;
	numtk::specfunc::bessel::Knu Knu;
	for (double x = 0.001; x < 2.0; x += 0.001)
		std::cout << x << ", " << Knu(1.0/3.0, x) << ", " << gamma/2.0*std::pow(2.0/x, 1.0/3.0) <<  std::endl;
	return 0;
}