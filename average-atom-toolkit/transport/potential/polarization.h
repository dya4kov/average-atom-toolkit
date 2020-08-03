#pragma once
#include <average-atom-toolkit/transport/potential/base.h>

namespace aatk {
namespace transport {
namespace potential {

class Polarization : public Base {
public:
	explicit Polarization(double r0, double alpha, double _eps = 1.e-6);
	double  operator()(double r);
	double  delta_eps(int l, double k) override;
private:
	double r0;
	double alpha;
};

}
}
}