#pragma once
#include <vector>
#include <cstddef>

#include <average-atom-toolkit/configuration.h>

namespace aatk {
namespace TF {
namespace shell {

class FreeEnergy {
public:
    FreeEnergy();

    double operator() (const double V, const double T);
    double DV         (const double V, const double T);
    double DT         (const double V, const double T);
    double D2V        (const double V, const double T);
    double DVT        (const double V, const double T);
    double D2T        (const double V, const double T);

    std::vector<double> operator() (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double> DV         (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double> DT         (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double> D2V        (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double> DVT        (const std::vector<double>& V, const std::vector<double>& T);
    std::vector<double> D2T        (const std::vector<double>& V, const std::vector<double>& T);

    double* operator() (const double* V, const double* T, const size_t vsize, const size_t tsize);
    double* DV         (const double* V, const double* T, const size_t vsize, const size_t tsize);
    double* DT         (const double* V, const double* T, const size_t vsize, const size_t tsize);
    double* D2V        (const double* V, const double* T, const size_t vsize, const size_t tsize);
    double* DVT        (const double* V, const double* T, const size_t vsize, const size_t tsize);
    double* D2T        (const double* V, const double* T, const size_t vsize, const size_t tsize);

    void setZ(const double Z);
    void setTolerance(const double eps);
    void setNmax(const int nmax);

#ifdef ENABLE_MULTITHREADING
    void setThreadsLimit(const std::size_t N);
#endif

private:

    double F   (const double V, const double T);
    double FDV (const double V, const double T);
    double FDT (const double V, const double T);
    double FD2V(const double V, const double T);
    double FDVT(const double V, const double T);
    double FD2T(const double V, const double T);
    
    double tolerance;
    double Z;
    int    nmax;

    static const std::size_t dnsh_interpolation_size = 3501;

#ifdef ENABLE_MULTITHREADING
    std::size_t threadsLimit;
#endif

};

}
}
}