#pragma once
#include <average-atom-toolkit/transport/potential/base.h>

namespace aatk {
namespace transport {
namespace potential {

class Debye : public Base {
public:
	explicit Debye(double Z, double _eps = 1.e-6);
	double  operator()(double r) override;
	double  delta_eps(int l, double k) override;
private:
	double Z;
};

}
}
}