#include <cmath>
#include <cstddef>
#include <average-atom-toolkit/transport/potential/base.h>

namespace aatk {
namespace transport {
namespace potential {

double* Base::operator()(double* r, double* V, std::size_t npoints) {
	auto& pot = *this;
	for (std::size_t i = 0; i < npoints; ++i) {
		V[i] = pot(r[i]);
	}
	return V;
}

}
}
}
