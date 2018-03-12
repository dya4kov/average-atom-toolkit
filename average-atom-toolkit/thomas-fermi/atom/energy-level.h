#pragma once
#include <vector>
#include <cstddef>
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
    // multithreaded call
    void prepareLevelsBelow(const int& nMax);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setVTZ(
        const double& V, 
        const double& T, 
        const double& Z
    );

    void setTolerance(const double& t);
    void setThreadsLimit(const std::size_t& N);

private:

    double V, T, Z;
    double tolerance;
    std::size_t threadsLimit;
    Action action;

    inline int iLevel(const int& n) const { return (n*(n - 1)/2); };

    void runLevel(const int& n, const int& l);
    void updateThreads(int& threads, int& current, int& last);
    ::std::function<void(const int& n, const int& l)> p_runLevel;

    ::std::vector<double>       eLevelBuffer;
    ::std::vector<bool>         eLevelReady;
    const ::std::vector<double> eLevelStart;

    bool bufferOk(const int& n);
    ::std::vector<int> needLevels(const int& n);
    
};

}
}