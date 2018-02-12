#pragma once
#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/atom/potential.h>
#include <average-atom-tools/thomas-fermi/atom/ODE/rotate-points.h>

namespace AATools {
namespace TF {

using ::numtools::ODE::Array;
using ::numtools::ODE::Solver;
using ::numtools::ODE::stepper::PD853;

using ODE::RHSRP1;
using ODE::RHSRP2;

class RotatePoints {
public:
    RotatePoints();
    RotatePoints(const RotatePoints& rps);
    RotatePoints& operator=(const RotatePoints& rps);

    double inner(const double& e, const double& l);
    double outer(const double& e, const double& l);

    Array<RHSRP1::dim> innerY(const double& e, const double& l);
    Array<RHSRP2::dim> outerY(const double& e, const double& l);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setTolerance(const double& t);

private:

    void setRP1();
    void setRP2();
    void setParam(const double& e, const double& l);

    double tolerance;

    double V1, T1, mu1;
    double VZ, TZ, muZ;

    double e, l;

    double rp1; bool rp1ready; 
    double rp2; bool rp2ready;

    Array<RHSRP1::dim> y1;
    Array<RHSRP2::dim> y2;

    RHSRP1 rhsRP1;
    RHSRP2 rhsRP2;

    Solver<PD853<RHSRP1>> solverRP1;
    Solver<PD853<RHSRP2>> solverRP2;

    Potential potential;

};

}
}