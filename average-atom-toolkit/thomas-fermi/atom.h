#pragma once
#include <vector>

#include <average-atom-toolkit/thomas-fermi/atom/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>
#include <average-atom-toolkit/thomas-fermi/atom/electron-density.h>

namespace aatk {
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

    void setVTZ(
        const double& V, 
        const double& T, 
        const double& Z
    );

    void setTolerance (const double& t);

private:

    double V, T, Z;
    double tolerance;

    RotatePoints RP;
};

}
}