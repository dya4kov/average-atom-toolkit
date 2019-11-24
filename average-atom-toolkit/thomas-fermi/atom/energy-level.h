#pragma once
#include <vector>
#include <utility>
#include <cstddef>
#include <functional>
#include <average-atom-toolkit/configuration.h>
#include <average-atom-toolkit/thomas-fermi/atom/action.h>

namespace aatk {
namespace TF {

class EnergyLevel {
public:
    EnergyLevel();
    EnergyLevel(const EnergyLevel& e);
    EnergyLevel& operator=(const EnergyLevel& e);
    // single thread call
    double operator()(const std::size_t n, const std::size_t l);
    // multithreaded call
    std::vector<double> operator[](const std::size_t n);
    // multithreaded call
    void prepareLevelsBelow(const std::size_t nMax);

    void setV(const double V);
    void setT(const double T);
    void setZ(const double Z);

    void setVTZ(
        const double V, 
        const double T, 
        const double Z
    );

    void setTolerance(const double t);

#ifdef ENABLE_MULTITHREADING
    void setThreadsLimit(const std::size_t N);
#endif

private:

    double V, T, Z;
    double tolerance;
    Action action;
    std::size_t nmax;

    void setNmax(const std::size_t nmax);
    void resetLevels();
    void runLevel(const std::size_t n, const std::size_t l);

    std::vector<std::pair<std::size_t,std::size_t>> needLevels();
    std::vector<std::pair<std::size_t,std::size_t>> needLevels(const std::size_t n);

    void runLevels(const std::vector<std::pair<std::size_t,std::size_t>>& tasks);
    std::function<void(const std::vector<std::pair<std::size_t,std::size_t>>& tasks)> p_runLevels;

#ifdef ENABLE_MULTITHREADING
    std::size_t threadsLimit;
#endif

    std::vector<std::vector<double>> eLevel; 
    std::vector<std::vector<double>> eReady;

    const std::vector<double> eLevelStart;

    
};

}
}