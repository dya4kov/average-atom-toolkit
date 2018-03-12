#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>

namespace numtk {
namespace specfunc {
namespace FDI {

HalfInc::HalfInc() {
    solver.setTolerance(1e-10, 1e-10);
}
double HalfInc::value(const double& x, const double& y) {
    Array<1> result; result.fill(0.0);    
    rhs.set_p(x);
    solver.integrate(rhs, result, 0.0, y);
    return result[0];
}

HalfInc2::HalfInc2() {
    solver.setTolerance(1e-10, 1e-10);
}

double HalfInc2::value(const double& x, const double& y) {
    double offset = 25.0;
    double result = 0.0;

    if (x - offset > 0.0) {
        result += 2.0/3.0*std::min(y, x - offset)*std::sqrt(std::min(y, x - offset));
    }
    if (y > x - offset) {
        Array<1> i; i.fill(0.0);
        rhs.set_p(x);
        solver.integrate(rhs, i, std::max(0.0, x - offset), std::min(y, x + offset));
        result += i[0];
    }
    if (y > x + offset) {
        result += std::exp(x)*(
          0.5*std::sqrt(M_PI)*(
              std::erf(std::sqrt(y)) - 
              std::erf(std::sqrt(std::max(0.0, x + offset)))
          )
          -
          (
              std::sqrt(y)*std::exp(-y) -
            std::sqrt(std::max(0.0, x + offset))*std::exp(-std::max(0.0, x + offset))
          )
        );
    }

    return result;
}

}
}
}