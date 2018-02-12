#pragma once
#include <numeric-tools/ODE/types.h>
#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/atom/potential.h>
#include <average-atom-tools/thomas-fermi/atom/rotate-points.h>
#include <average-atom-tools/thomas-fermi/atom/ODE/action.h>

namespace AATools {
namespace TF {

using ::numtools::ODE::Solver;
using ::numtools::ODE::stepper::PD853;

using ODE::RHSAction;

class Action {
public:
    Action();
    Action(const Action& a);
    Action& operator=(const Action& a);

    double operator()(const double& e, const double& l);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setTolerance(const double& t);
private:

    void setAction();
    void setParam(const double& e, const double& l);

    double tolerance;

    double V1, T1, mu1;
    double VZ, TZ, muZ;

    double e, l;

    double action; 
    bool   ready;

    RHSAction rhs;
    Solver<PD853<RHSAction>> solver;

    RotatePoints RP;
    Potential potential;
};

}
}