#pragma once
#include <average-atom-toolkit/transport/potential/base.h>

namespace aatk {
namespace transport {
namespace potential {

class SquareWell : public Base {
public:
	explicit SquareWell(double width, double depth, double _eps = 1.e-6);
	double  operator()(double r);
	double  delta_eps(int l, double k) override;
private:
	double width;
	double depth;
};

}
}
}