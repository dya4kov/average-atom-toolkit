#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>

namespace numtk {
namespace specfunc {
namespace FDI {

ThreeHalfInc::ThreeHalfInc() {
    solver.setTolerance(1e-10, 1e-10);
}
double ThreeHalfInc::value(const double& x, const double& y) {
    Array<1> result; result.fill(0.0);    
    rhs.set_p(x);
    solver.integrate(rhs, result, 0.0, y);
    return result[0];
}

}
}
}