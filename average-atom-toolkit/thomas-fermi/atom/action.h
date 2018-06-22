#pragma once
#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>

namespace aatk {
namespace TF {

class Action {
public:
    Action();
    Action(const Action& a);
    Action& operator=(const Action& a);

    double operator()(const double& energy, const double& lambda);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setVTZ(
        const double& V, 
        const double& T, 
        const double& Z
    );

    void setTolerance(const double& t);
private:

    void setAction();
    void setParam(const double& energy, const double& lambda);

    double tolerance;
    double V, T, Z, mu;
    double energy, lambda;

    double action; 
    bool   ready;

    RotatePoints RP;
};

}
}