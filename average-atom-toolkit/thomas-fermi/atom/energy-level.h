#pragma once
#include <vector>
#include <functional>
#include <average-atom-toolkit/thomas-fermi/atom/action.h>

namespace aatk {
namespace TF {

class EnergyLevel {
public:
    EnergyLevel();
    EnergyLevel(const EnergyLevel& e);
    EnergyLevel& operator=(const EnergyLevel& e);
    // single thread call
    double operator()(const int& n, const int& l);
    // multithreaded call
    std::vector<double> operator[](const int& n);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setTolerance(const double& t);

private:

    double V, T, Z;
    double tolerance;
    Action action;

    inline int iLevel(const int& n) const { return (n*(n - 1)/2); };

    void runLevel(const int& n, const int& l);
    ::std::function<void(const int& n, const int& l)> p_runLevel;

    ::std::vector<double>       eLevelBuffer;
    ::std::vector<bool>         eLevelReady;
    const ::std::vector<double> eLevelStart;

    void checkBufSize(const int& n);
    ::std::vector<int> needLevels(const int& n);
    
};

}
}