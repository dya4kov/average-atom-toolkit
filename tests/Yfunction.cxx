#include <iostream>
#include <numeric-tools/specfunc/fermi-dirac/Yfunction.h>

using namespace numtools::specfunc;

int main() {

	Yfunction Y;

	int nx = 100;

    std::vector<double> x(nx);

    double xmin = -1.0; 
    double xmax =  1.0;

    for (int ix = 0; ix < nx; ++ix)  {
        x[ix] = xmin + ix*(xmax - xmin)/(nx - 1);
    }

	auto y = Y(x);

	for (auto&& yi : y) std::cout << yi << std::endl;

	auto dy = Y.derivative(x);

	for (auto&& yi : dy) std::cout << yi << std::endl;

	return 0;
}