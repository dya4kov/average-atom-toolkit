#include <cmath>
#include <cstddef>
#include <average-atom-toolkit/transport/cross-section.h>

namespace aatk {
namespace transport {

CrossSection::CrossSection(
	potential::Base& pot, 
	int             _nterms, 
	double          _rmax) : 
    delta(pot, _rmax),  
    nterms(_nterms) {}

double CrossSection::operator()(double k) {
	double sigma;
    double deltaPrev, deltaNext, diff;
    deltaPrev = delta(0, k);
    sigma = 0.;
    for(int l = 0; l < nterms; ++l) {
        deltaNext = delta(l + 1, k);
        diff = std::sin(deltaPrev - deltaNext);
        sigma += (l + 1) * diff * diff;
        deltaPrev = deltaNext;
    }
    return 4. * M_PI * sigma / k / k;
}

double* CrossSection::operator()(double* k, double* sigma, std::size_t npoints) {
	auto& cs = *this;
	for (std::size_t i = 0; i < npoints; ++i) {
		sigma[i] = cs(k[i]);
	}
	return sigma;
}

}
}
