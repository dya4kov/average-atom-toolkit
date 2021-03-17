#pragma once
#include <vector>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>

#include <average-atom-toolkit/util/specfunc/fermi-dirac/complete.h>

namespace aatk {
namespace util {
namespace specfunc {

using FD::DMHalf;
using FD::MHalf;
using FD::Half;

class Yfunction {
public:
	Yfunction(double tolerance = 1e-7);
    ~Yfunction();

    double operator()(const double x);
    double derivative(const double x);
    
    std::vector<double>& operator()(const std::vector<double>& x);
    std::vector<double>& derivative(const std::vector<double>& x);

    double* operator()(const double* x, const size_t& n);
    double* derivative(const double* x, const size_t& n);

private:

    double integral(const double x);
	double tolerance;

    FermiDirac<DMHalf> FDdmhalf;
    FermiDirac<MHalf>  FDmhalf;
    FermiDirac<Half>   FDhalf;

    gsl_odeiv2_step    *stepper;
    gsl_odeiv2_control *control;
    gsl_odeiv2_evolve   *evolve;
};

}
}
}