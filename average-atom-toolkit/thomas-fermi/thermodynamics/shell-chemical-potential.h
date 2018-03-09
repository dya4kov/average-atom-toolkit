#pragma once
#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/chemical-potential.h>

namespace aatk {
namespace TF {
namespace shell {

class ChemicalPotential {
public:
    ChemicalPotential();
    ChemicalPotential(const ChemicalPotential& rps);
    ChemicalPotential& operator=(const ChemicalPotential& rps);

    double operator()(const double& V, const double& T);

    void setZ(const double& Z);
    void setThreadsLimit(const std::size_t& N);
    void setTolerance(const double& t);

private:

    double tolerance;
    double Z;
    std::size_t threadsLimit;

    void M(const double& V, const double& T, double& result, bool& finished);

    double* evaluate(
        ::std::function<void(const double&, const double&, double&, bool&)> func,
          const double* V, 
          const double* T, 
          const std::size_t& vsize, 
          const std::size_t& tsize
    );

    ::aatk::TF::ChemicalPotential mu;
    ::aatk::TF::ElectronStates     N;

};

}
}
}