#pragma once
#include <vector>

#include <average-atom-tools/thomas-fermi/atom/potential.h>
#include <average-atom-tools/thomas-fermi/atom/rotate-points.h>
#include <average-atom-tools/thomas-fermi/atom/energy-level.h>
#include <average-atom-tools/thomas-fermi/atom/electron-states.h>
#include <average-atom-tools/thomas-fermi/atom/electron-density.h>

namespace AATools {
namespace TF {

class Atom {
public:
    Atom();
    Atom(const Atom& a);
    Atom& operator=(const Atom& a);

    Potential          xU;
    EnergyLevel         e;
    ElectronStates      N;
    ElectronDensity eDens;

    double rpInner(const int& n, const int& l);
    double rpOuter(const int& n, const int& l);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

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