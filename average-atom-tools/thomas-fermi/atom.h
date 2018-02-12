#pragma once
#include <vector>

#include <average-atom-tools/thomas-fermi/atom/potential.h>
#include <average-atom-tools/thomas-fermi/atom/rotate-points.h>
#include <average-atom-tools/thomas-fermi/atom/energy-level.h>

namespace AATools {
namespace TF {

class Atom {
public:
    Atom();
    Atom(const Atom& a);
    Atom& operator=(const Atom& a);

    Potential  xU;
    EnergyLevel e;

    double               eDens (const double& x);
    double*              eDens (double* x, size_t n);
    std::vector<double>& eDens (const std::vector<double>& x);

    double rpInner(const int& n, const int& l);
    double rpOuter(const int& n, const int& l);

    void setV(const double& _V);
    void setT(const double& _T);
    void setZ(const double& _Z);

    void setTolerance (const double& t);

private:

    double V1,  VZ;
    double T1,  TZ;
    double mu1, muZ;
    double tolerance;
    
    RotatePoints RP;
};

}
}