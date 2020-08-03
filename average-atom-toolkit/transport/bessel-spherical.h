#pragma once
#include <numeric-toolkit/specfunc/bessel/Jnu.h>
#include <numeric-toolkit/specfunc/bessel/Ynu.h>

namespace aatk {
namespace transport {

struct Jl {
	double operator()(const int l, const double x) {
		return std::sqrt(0.5 * M_PI * x) * gsl_sf_bessel_Jnu(l + 0.5, x);
	}
};

struct Nl {
	double operator()(const int l, const double x) {
		return std::sqrt(0.5 * M_PI * x) * gsl_sf_bessel_Ynu(l + 0.5, x);
	}
};

}
}