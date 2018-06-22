#pragma once

namespace aatk {
namespace TF {

class RotatePoints {
public:
    RotatePoints();
    RotatePoints(const RotatePoints& rps);
    RotatePoints& operator=(const RotatePoints& rps);

    double inner(const double& energy, const double& lambda);
    double outer(const double& energy, const double& lambda);

    double* innerY(const double& energy, const double& lambda);
    double* outerY(const double& energy, const double& lambda);

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

    void setRPinner();
    void setRPouter();
    void setParam(const double& energy, const double& lambda);

    double tolerance;
    double V, T, Z, mu;
    double energy, lambda;

    double rpInner; bool rpIready; 
    double rpOuter; bool rpOready;
    double yOuter[2];
    double yInner[2];

};

}
}