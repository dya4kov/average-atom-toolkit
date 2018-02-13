#pragma once
#include <vector>
#include <functional>

namespace AATools {
namespace TF {

class FreeEnergy {
public:
    FreeEnergy();

    double operator() (const double& V, const double& T);
    double DV         (const double& V, const double& T);
    double DT         (const double& V, const double& T);
    double D2V        (const double& V, const double& T);
    double DVT        (const double& V, const double& T);
    double D2T        (const double& V, const double& T);

    std::vector<double>& operator() (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double>& DV         (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double>& DT         (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double>& D2V        (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double>& DVT        (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double>& D2T        (const std::vector<double>& V, const std::vector<double>& T);

    double* operator() (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    double* DV         (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    double* DT         (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    double* D2V        (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    double* DVT        (const double* V, const double* T, const size_t& vsize, const size_t& tsize);
    double* D2T        (const double* V, const double* T, const size_t& vsize, const size_t& tsize);

    void setZ(const double& Z);
    void setTolerance(const double& eps);

private:

    void F   (const double& V, const double& T, double& result, bool& finished);
    void FDV (const double& V, const double& T, double& result, bool& finished);
    void FDT (const double& V, const double& T, double& result, bool& finished);
    void FD2V(const double& V, const double& T, double& result, bool& finished);
    void FDVT(const double& V, const double& T, double& result, bool& finished);
    void FD2T(const double& V, const double& T, double& result, bool& finished);

	double* evaluate(::std::function<void(const double&, const double&, double&, bool&)> func,
		           const double* V, const double* T, const size_t& vsize, const size_t& tsize);

    double tolerance;
    double Z;
    const double E0;

};

}
}