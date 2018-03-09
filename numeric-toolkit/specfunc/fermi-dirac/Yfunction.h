#pragma once
#include <vector>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/Yfunction/ODE.h>

namespace numtk {
namespace specfunc {
    
using ::numtk::ODE::Solver;
using ::numtk::ODE::stepper::PD853;

using ODE::RHSY;

using FD::DMHalf;
using FD::MHalf;
using FD::Half;

class Yfunction {
public:
	Yfunction();

    double operator()(const double& x);
    double derivative(const double& x);
    
    std::vector<double>& operator()(const std::vector<double>& x);
    std::vector<double>& derivative(const std::vector<double>& x);

    double* operator()(const double* x, const size_t& n);
    double* derivative(const double* x, const size_t& n);

    void setTolerance(const double& eps);

private:

    double integral(const double& x);

	double tolerance;

    FermiDirac<DMHalf> FDdmhalf;
    FermiDirac<MHalf>  FDmhalf;
    FermiDirac<Half>   FDhalf;

    RHSY rhs;
    Solver<PD853<RHSY>> solver;

};
}
}