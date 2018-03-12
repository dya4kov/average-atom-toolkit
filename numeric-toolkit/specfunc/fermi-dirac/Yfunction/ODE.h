#pragma once
#include <cmath>
#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

namespace numtk {
namespace specfunc {
namespace ODE {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

using ::numtk::specfunc::FermiDirac;
using ::numtk::specfunc::FD::MHalf;

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