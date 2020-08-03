#include <cmath>
#include <average-atom-toolkit/transport/potential/polarization.h>

namespace aatk {
namespace transport {
namespace potential {

Polarization::Polarization(double _r0, double _alpha, double _eps) :
                               r0(_r0), alpha(_alpha),  Base(_eps)
{}

double Polarization::delta_eps(int l, double k) {
	double df = doubleFact(2 * l + 1); df *= df;
    return 0.5 * alpha / std::pow(r0, 4.) * std::pow(k * eps, 2. * l + 1) * 
           eps * eps / (2. * l + 3.) / df;
}

double Polarization::operator()(double r) {
	double denom = r * r + r0 * r0; 
    return -0.5 * alpha / denom / denom;
}

}
}
}
