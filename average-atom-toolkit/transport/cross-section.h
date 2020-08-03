#pragma once
#include <vector>
#include <cstddef>
#include <average-atom-toolkit/transport/potential/base.h>
#include <average-atom-toolkit/transport/delta.h>

namespace aatk {
namespace transport {

class CrossSection {
public:
	CrossSection(potential::Base& pot, int nterms = 30, double rmax = 10000.);
	double  operator()(double k);
	double* operator()(double* k, double* sigma, std::size_t npoints);
private:
	Delta delta;
	int  nterms;
};

}
}
