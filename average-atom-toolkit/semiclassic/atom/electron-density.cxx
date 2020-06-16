#include <cmath>
#include <average-atom-toolkit/semiclassic/atom.h>

namespace aatk {
namespace semiclassic {

std::vector<double> Atom::electronDensity(const std::vector<double>& x) {
	std::vector<double> result(x.size(), 0.0);
	electronDensity(x.data(), result.data(), x.size());
	return result;
}
double Atom::electronDensity(double x) {
	double result;
	electronDensity(&x, &result, 1);
	return result;
}
void Atom::electronDensity(const double* x, double* y, std::size_t n) {
	auto& density = *densityInterpolation;
	for (std::size_t i = 0; i < n; ++i) {
		y[i] = density(std::sqrt(x[i]));
	}
}

}
}