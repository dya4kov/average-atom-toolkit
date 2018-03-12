#pragma once
#include <vector>
#include <functional>

namespace aatk {
namespace TF {

class ChemicalPotential {
public:

    ChemicalPotential();
    ChemicalPotential(const ChemicalPotential& mu);
    ChemicalPotential& operator=(const ChemicalPotential& mu);

    double operator()(const double& VZ, const double& TZ);
    std::vector<double> operator() (const std::vector<double>& VZ, const std::vector<double>& TZ);
    double* operator()(const double* VZ, const double* TZ, const std::size_t& vsize, const std::size_t& tsize);

    void setZ(const double& Z);
    void setThreadsLimit(const std::size_t& N);
    void setTolerance(const double& eps);

private:

    void   mu(const double& V, const double& T, double& result, bool& finished);

    double mu1(const double& V1, const double& T1, const double& tol);

    double mu1_approx(const double& lgV1  /*    T = 0     */);
    double mu1_approx(const double& lgV1, const double& lgT1);

    void updateThreads(
        std::size_t& threads, 
        std::size_t& current, 
        std::size_t& last, 
        bool* finished
    );

    std::vector<double> evaluate(
        ::std::function<void(const double&, const double&, double&, bool&)> func,
          const double* V, 
          const double* T,
          const std::size_t& vsize, 
          const std::size_t& tsize
    );

    // double localTolerance;
    // double localZ;
    double Z;
    double tolerance;
    std::size_t threadsLimit;

    // 
    // static double globalTolerance;
    // static double globalZ;
    

    // current buffer
    // static std::mutex       mutexBuf;
    // static std::vector<double>  vBuf;
    // static std::vector<double>  tBuf;
    // static std::vector<double> muBuf;

    // bool validBuf();

    // int vIndex(const double& V);
    // int tIndex(const double& T);
    // void prepareBuf(
    //     const double* V, 
    //     const double* T, 
    //     const std::size_t& vsize, 
    //     const std::size_t& tsize
    // );
    // std::vector<double>  sorted(
    //     const double*     array, 
    //     const std::size_t& size
    // );

    // basic buffer
    static std::vector<double>   table;
    static const int    vSize,   tSize;
    static const double lgV0,    lgT0;
    static const double lgVstep, lgTstep;
    static const double bestTolerance;
};

}
}
