#pragma once
#include <cmath>

#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

namespace numtools {
namespace specfunc {
/**
* Incomplete Fermi-Dirac functions
*/
template<typename Func>
struct FermiDiracInc {
	FermiDiracInc() : f() {}
	double operator() (const double& x, const double& y) {
		return f.value(x, y);
	}
private:
	Func f;
};

namespace FDI {

using ::numtools::ODE::Array;
using ::numtools::ODE::Dimension;

using ::numtools::ODE::Solver;
using ::numtools::ODE::stepper::PD853;

class HalfInc {
public:
	HalfInc();
	double value(const double& x, const double& y);
private:	
	struct RHS {
		static const Dimension dim = 1;
		RHS() {}
		void set_p(const double& _p) { p = _p; }
		void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) {
			dydx[0] = std::sqrt(x)/(1.0 + std::exp(x - p));
		}
	private:
		double p;
	};

	RHS rhs; Solver<PD853<RHS>> solver;
};

class HalfInc2 {
public:
	HalfInc2();
	double value(const double& x, const double& y);
private:	
	struct RHS {
		static const Dimension dim = 1;
		RHS() {}
		void set_p(const double& _p) { p = _p; }
		void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) {
			dydx[0] = std::sqrt(x)/(1.0 + std::exp(x - p));
		}
	private:
		double p;
	};

	RHS rhs; Solver<PD853<RHS>> solver;
};
}
}
}