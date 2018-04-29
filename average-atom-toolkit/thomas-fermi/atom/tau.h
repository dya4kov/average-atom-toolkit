#pragma once
#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>

namespace aatk {
namespace TF {

class Tau {
public:
    Tau();
    Tau(const Tau& a);
    Tau& operator=(const Tau& a);

    double operator()(const double& e, const double& l);

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

    void setTau();
    void setParam(const double& e, const double& l);

    double tolerance;
    double V, T, Z, mu;
    double e, l;

    double tau; 
    bool   ready;

    RotatePoints RP;
};

}
}