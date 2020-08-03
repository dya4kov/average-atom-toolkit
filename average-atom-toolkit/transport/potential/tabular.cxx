#include <algorithm>
#include <numeric>
#include <cmath>

#include <average-atom-toolkit/transport/potential/tabular.h>

namespace aatk {
namespace transport {
namespace potential {

Tabular::Tabular(double* _r, double* _rV, std::size_t npoints, double _eps) : Base(_eps)
{
	std::vector<double> sqr(npoints);
    std::vector<double> rV(npoints);
	bool needSort = false;
	for (std::size_t i = 0; i < npoints - 1; ++i) {
		needSort = needSort || _r[i] > _r[i + 1];
	}
	if (needSort) {
		std::vector<std::size_t> idx(npoints);
    	std::iota(idx.begin(), idx.end(), 0);
    	std::sort(idx.begin(), idx.end(),
    	   [&_r](std::size_t i1, std::size_t i2) {
    	        return _r[i1] < _r[i2];
    	   }
    	);
    	std::size_t k = 0;
    	for (auto i : idx) {
			sqr[k] = std::sqrt(_r[i]);
			rV[k]  = _rV[i];
			++k;
    	}
	}
	else {
		for (std::size_t i = 0; i < npoints; ++i) {
			sqr[i] = std::sqrt(_r[i]);
			rV[i]  = _rV[i];
    	}
	}
	rmin = sqr.front()*sqr.front();
	rmax = sqr.back()*sqr.back();
	rVmin = rV.front();
	rVmax = rV.back();
    interpolation = new Spline(sqr, rV);
}
Tabular::~Tabular() {
	delete interpolation;
}

double Tabular::operator()(double r) {
	auto& potential = *interpolation;
	return r < rmax ? potential(std::sqrt(r))/r : 0.;
}

double Tabular::delta_eps(int l, double k) {
	double df = doubleFact(2 * l + 1);
	df *= df;
    return -rVmin*std::pow(k * eps, 2.*l + 1.)*eps / (2.*l + 2) / df;
}

}
}
}
