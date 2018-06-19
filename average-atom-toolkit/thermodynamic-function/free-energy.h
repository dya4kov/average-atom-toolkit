#pragma once
#include <cstddef>
#include <string>
#include <vector>

#include <average-atom-toolkit/thermodynamic-function/thermodynamic-function.h>

namespace aatk {
namespace tfunc {

class FreeEnergy : ThermodynamicFunction {
public:

    virtual double operator() (const double& V, const double& T);
    virtual double DV         (const double& V, const double& T);
    virtual double DT         (const double& V, const double& T);
    virtual double D2V        (const double& V, const double& T);
    virtual double DVT        (const double& V, const double& T);
    virtual double D2T        (const double& V, const double& T);

    virtual std::vector<double> operator() (const std::vector<double>& V, const std::vector<double>& T);
    virtual std::vector<double> DV         (const std::vector<double>& V, const std::vector<double>& T);
    virtual std::vector<double> DT         (const std::vector<double>& V, const std::vector<double>& T);
    virtual std::vector<double> D2V        (const std::vector<double>& V, const std::vector<double>& T);
    virtual std::vector<double> DVT        (const std::vector<double>& V, const std::vector<double>& T);
    virtual std::vector<double> D2T        (const std::vector<double>& V, const std::vector<double>& T);

    virtual double* operator() (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    virtual double* DV         (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    virtual double* DT         (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    virtual double* D2V        (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    virtual double* DVT        (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    virtual double* D2T        (const double* V, const double* T, const size_t& vsize, const size_t& tsize);

    std::string name() { return std::string("FreeEnergy"); }

};

}
}