#pragma once
#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>

namespace aatk {
namespace TF {

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
    double V, T, Z, mu;
    double e, l;

    double action; 
    bool   ready;

    RotatePoints RP;
};

}
}