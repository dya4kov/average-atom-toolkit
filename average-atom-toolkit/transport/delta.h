#pragma once
#include <vector>
#include <cstddef>
#include <average-atom-toolkit/transport/potential/base.h>

namespace aatk {
namespace transport {

class Delta {
public:
	Delta(potential::Base& pot, double rmax = 10000.);
	double operator()(int l, double k);
private:
	potential::Base& potential;
	double           rmax;
};

}
}
