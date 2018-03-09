#pragma once

namespace aatk {
namespace TF {

class RotatePoints {
public:
    RotatePoints();
    RotatePoints(const RotatePoints& rps);
    RotatePoints& operator=(const RotatePoints& rps);

    double inner(const double& e, const double& l);
    double outer(const double& e, const double& l);

    double* outerY(const double& e, const double& l);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setTolerance(const double& t);

private:

    void setRPinner();
    void setRPouter();
    void setParam(const double& e, const double& l);

    double tolerance;
    double V, T, Z, mu;
    double e, l;

    double rpInner; bool rpIready; 
    double rpOuter; bool rpOready;
    double yOuter[2];

};

}
}