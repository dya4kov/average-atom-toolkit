#include <cmath>
#include <iostream>
#include <average-atom-toolkit/transport/potential/debye.h>

namespace aatk {
namespace transport {
namespace potential {

Debye::Debye(double _Z, double _eps) : 
                  Z(_Z),  Base(_eps) {}

double Debye::delta_eps(int l, double k) {
	double df = doubleFact(2 * l + 1); df *= df;
    return Z * std::pow(k * eps, 2. * l + 1.) * eps / (2. * l + 2.) / df;
}

double Debye::operator()(double r) {
	return -Z / r * std::exp(-r);
}

}
}
}
