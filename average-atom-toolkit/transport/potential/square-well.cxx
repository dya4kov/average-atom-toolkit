#include <cmath>
#include <average-atom-toolkit/transport/potential/square-well.h>

namespace aatk {
namespace transport {
namespace potential {

SquareWell::SquareWell(double _width, double _depth, double _eps) :
                        width(_width), depth(_depth),  Base(_eps)
{}

double SquareWell::delta_eps(int l, double k) {
	double df = doubleFact(2 * l + 1); df *= df;
    return depth * std::pow(k * eps, 2. * l + 1.) * eps * eps / (2. * l + 3.) / df;
}

double SquareWell::operator()(double r) {
	return r < width ? -depth : 0.;
}

}
}
}
