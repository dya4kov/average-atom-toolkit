#pragma once
#include <numeric-toolkit/interpolation/spline.h>
#include <average-atom-toolkit/transport/potential/base.h>

namespace aatk {
namespace transport {
namespace potential {

using ::numtk::interpolation::Spline;

class Tabular : public Base {
public:
	Tabular(double* r, double* V, std::size_t npoints, double _eps = 1.e-6);
	~Tabular();

	double operator()(double r) override;
	double delta_eps(int l, double k) override;

private:
	Spline* interpolation;
	double  rmin, rVmin;
	double  rmax, rVmax;
};

}
}
}