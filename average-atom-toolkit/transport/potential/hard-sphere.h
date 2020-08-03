#pragma once
#include <average-atom-toolkit/transport/potential/base.h>

namespace aatk {
namespace transport {
namespace potential {

class HardSphere : public Base {
public:
	explicit HardSphere(double _eps = 1.e-6);
	double  operator()(double r);
	double  delta_eps(int l, double k) override;
};

}
}
}