#pragma once
#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

namespace numtk {
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

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::ODE::Solver;
using ::numtk::ODE::stepper::PD853;

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