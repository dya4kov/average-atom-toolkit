#pragma once
#include <cmath>
#include <numeric-tools/ODE/types.h>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

namespace numtools {
namespace specfunc {
namespace ODE {

using ::numtools::ODE::Array;
using ::numtools::ODE::Dimension;

using ::numtools::specfunc::FermiDirac;
using ::numtools::specfunc::FD::MHalf;

struct RHSY {
    public:
    static const Dimension dim = 1;
    static const Dimension result = 0;
    RHSY() {}
    void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) {
        double f = FDmhalf(x);
        dydx[0] = f*f;
    }
    private:
    FermiDirac<MHalf> FDmhalf;
};    

}
}
}