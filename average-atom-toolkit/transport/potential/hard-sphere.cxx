#include <cmath>
#include <average-atom-toolkit/transport/potential/hard-sphere.h>

namespace aatk {
namespace transport {
namespace potential {

HardSphere::HardSphere(double _eps) : 
                         Base(_eps) {}

double HardSphere::delta_eps(int l, double k) {
	return std::atan(Jl(l, k * eps) / Nl(l, k * eps));
}

double HardSphere::operator()(double r) {
	return r < eps ? 1.e200 : 0.;
}

}
}
}
