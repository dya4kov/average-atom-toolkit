#include <cmath>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/action.h>

#ifdef ENABLE_MULTITHREADING
#include <thread>
#endif

using namespace   aatk::TF;
using ::std::placeholders::_1;

EnergyLevel::EnergyLevel() :
    eLevelStart({-1e+3, -1e+2, -1e+1, -1.0, 0.0, 1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5}),
    V(1.0), T(1.0), Z(1.0), nmax(0),
#ifdef ENABLE_MULTITHREADING
    threadsLimit(std::max(4u, std::thread::hardware_concurrency())),
#endif
    tolerance(1e-6)
{
    p_runLevels = std::bind(&EnergyLevel::runLevels, this, _1);
}

EnergyLevel::EnergyLevel(const EnergyLevel& e) {
    V = e.V; T = e.T; Z = e.Z;

    action      = e.action;
    tolerance   = e.tolerance;
    nmax        = e.nmax;
    eLevel      = e.eLevel;
    eReady      = e.eReady;
    p_runLevels = std::bind(&EnergyLevel::runLevels, this, _1);

    return;
}

EnergyLevel& EnergyLevel::operator=(const EnergyLevel& e) {
    V = e.V; T = e.T; Z = e.Z;

    action      = e.action;
    tolerance   = e.tolerance;
    nmax        = e.nmax;
    eLevel      = e.eLevel;
    eReady      = e.eReady;
    p_runLevels = std::bind(&EnergyLevel::runLevels, this, _1);

    return *this;
}

void EnergyLevel::setV(const double _V) {
    V = _V;
    action.setV(V);
    resetLevels();
}

void EnergyLevel::setT(const double _T) {
    T = _T;
    action.setT(T);
    resetLevels();
}

void EnergyLevel::setZ(const double _Z) {
    Z = _Z;
    action.setZ(Z);
    resetLevels();
}

void EnergyLevel::setVTZ(
    const double _V,
    const double _T,
    const double _Z
) {
    V = _V; T = _T; Z = _Z;
    action.setVTZ(V,T,Z);
    resetLevels();
}

void EnergyLevel::setTolerance(const double t) {
    tolerance = t;
    action.setTolerance(t);
    resetLevels();
}

void EnergyLevel::resetLevels() {
    for (std::size_t n = 1; n <= nmax; ++n) {
        for (std::size_t l = 0; l < n; ++l) {
            eReady[n][l] = false;
        }
    }
}

#ifdef ENABLE_MULTITHREADING

void EnergyLevel::setThreadsLimit(const std::size_t Nthreads) {
    threadsLimit = Nthreads; 
}

#endif

void EnergyLevel::setNmax(const std::size_t _nmax) {
    if (nmax == _nmax) return;
    eLevel.resize(_nmax + 1); eReady.resize(_nmax + 1);
    eLevel[0].resize(0);  eReady[0].resize(0);
    for (std::size_t n = 1; n <= _nmax; ++n) {
        eLevel[n].resize(n);
        if (n <= nmax) eReady[n].resize(n);
        else eReady[n].resize(n, false);
    }
    nmax = _nmax;
}

void EnergyLevel::prepareLevelsBelow(const std::size_t _nmax) {
    setNmax(_nmax);
    auto levels = needLevels();

    if (levels.size() > 0) {

#ifdef ENABLE_MULTITHREADING
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> 
        tasks(threadsLimit);

    std::size_t ithread = 0;
    for (auto&& lvl : levels) {
        tasks[ithread].push_back(lvl);
        ithread = (ithread + 1) % threadsLimit;
    }

    std::vector<std::thread> threads;
    for (auto&& task : tasks) {
        threads.push_back(std::thread(p_runLevels, task));
    }
    for (auto&& thread : threads) {
        thread.join();
    }

#else
    for (auto&& lvl : levels) {
        auto n = lvl.first; auto l = lvl.second; runLevel(n, l);
    }
#endif

    } // if (levels.size() > 0)

    return;
}

std::vector<double> EnergyLevel::operator[](const std::size_t n) {
    if (n > nmax) setNmax(n);
    auto levels = needLevels(n);

    if (levels.size() > 0) {

#ifdef ENABLE_MULTITHREADING
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> 
        tasks(threadsLimit);

    std::size_t ithread = 0;
    for (auto&& lvl : levels) {
        tasks[ithread].push_back(lvl);
        ithread = (ithread + 1) % threadsLimit;
    }

    std::vector<std::thread> threads;
    for (auto&& task : tasks) {
        threads.push_back(std::thread(p_runLevels, task));
    }
    for (auto&& thread : threads) {
        thread.join();
    }
#else

    for (auto&& lvl : levels) {
        auto n = lvl.first; auto l = lvl.second; runLevel(n, l);
    }

#endif

    } // if (levels.size() > 0)

    std::vector<double> result(n);
    for (int l = 0; l < n; ++l)
        result[l] = std::pow(Z, 4.0/3.0)*eLevel[n][l];

    return result;
}

double EnergyLevel::operator()(const std::size_t n, const std::size_t l) {
    if (n > nmax) setNmax(n);
    if (!eReady[n][l]) runLevel(n, l);
    return std::pow(Z, 4.0/3.0)*eLevel[n][l];
}

std::vector<std::pair<std::size_t,std::size_t>> EnergyLevel::needLevels() {
    std::vector<std::pair<std::size_t,std::size_t>> need(nmax*(nmax + 1)/2);
    need.resize(0);
    for (std::size_t n = 1; n <= nmax; ++n) {
        for (std::size_t l = 0; l < n; ++l) {
            if (!eReady[n][l])
                need.push_back(std::pair<std::size_t,std::size_t>(n,l));
        }
    }
    return need;
}

std::vector<std::pair<std::size_t,std::size_t>> EnergyLevel::needLevels(const std::size_t n) {
    std::vector<std::pair<std::size_t,std::size_t>> need(n);
    need.resize(0);
    for (std::size_t l = 0; l < n; ++l) {
        if (!eReady[n][l])
            need.push_back(std::pair<std::size_t,std::size_t>(n,l));
    }
    return need;
}

// void EnergyLevel::runLevel(const int& n, const int& l) {
//     
//     Action act   = action;    
//     double exact = M_PI*(n - l - 0.5);
//     double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);
//     double lambda = l + 0.5;
//     double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
// 
//     double eLeft;
//     double eRight;
// 
//     int ie = 0;
//     while (act(eLevelStart[ie], lArg) < exact && ie < eLevelStart.size()) ++ie;
// 
//     eLeft  = eLevelStart[ie - 1];
//     eRight = eLevelStart[ie];
// 
//     if (eLevelReady[0]) eLeft = eLevelBuffer[0];
// 
//     int nStep = 0;
//     double err = std::abs(eLeft - eRight) / std::abs(eLeft + eRight);
//     while (err > 0.1*tolerance && nStep < 100) {
//         double eArg = 0.5*(eRight + eLeft);
//         if (act(eArg, lArg) - exact > 0.0) {
//             eRight -= 0.5*(eRight - eLeft);
//         }
//         else {
//             eLeft += 0.5*(eRight - eLeft);
//         }
//         err = std::abs(eLeft - eRight) / std::abs(eLeft + eRight);
//         ++nStep;
//     }
// 
//     // if (nStep == 100 ) {
//     //     std::cout << "Warning: energy level convergence failed" << std::endl;
//     // }
// 
//     //if (std::abs(eMin - eRight) < 10*tolerance) {
//     //    std::cout << "Warning: energy range is small" << std::endl;
//     //    std::cout << "Eleft = " << eMin << ", Eright = " << eMax << std::endl;
//     //}
// 
//     //if (std::abs(eMax - eLeft) < 10*tolerance) {
//     //    std::cout << "Warning: energy range is small" << std::endl;
//     //    std::cout << "Eleft = " << eMin << ", Eright = " << eMax << std::endl;
//     //}
// 
//     eLevelBuffer[iLevel(n) + l] = 0.5*(eLeft + eRight);
//     eLevelReady[iLevel(n) + l]  = true;
// }

void EnergyLevel::runLevels(const std::vector<std::pair<std::size_t,std::size_t>>& levels) {
    for (auto&& lvl : levels) {
        auto n = lvl.first;
        auto l = lvl.second;
        runLevel(n, l);
    }
}

void EnergyLevel::runLevel(const std::size_t n, const std::size_t l) {
    Action act   = action;
    double exact = M_PI*(n - l - 0.5);
    double lambda = l + 0.5;
    double Z43 = std::pow(Z, 4.0/3.0);
    int ie = 0;
    while (act(Z43*eLevelStart[ie], lambda) < exact && ie < eLevelStart.size()) ++ie;

    double eLeft  = eLevelStart[ie - 1];
    double eRight = eLevelStart[ie];

    // if (eLevelReady[0]) eLeft = eLevelBuffer[0];

    int nStep = 0;
    double dact = 1.e+5;
    double dactOld = 0.0;
    while (std::abs(dact - dactOld) > 0.1 || std::abs(dact - dactOld) < tolerance) {
        dactOld = dact;
        double eArg = 0.5*(eRight + eLeft);
        dact = act(Z43*eArg, lambda) - exact;
        if (dact > 0.0) {
            eRight -= 0.5*(eRight - eLeft);
        }
        else {
            eLeft += 0.5*(eRight - eLeft);
        }
//        std::cout << "nStep = " << nStep << ", dact = " << dact << ", dactOld = " << dactOld << std::endl;
        ++nStep;
    }

    double ePrev = eLeft;
    double dactPrev = act(Z43*ePrev, lambda) - exact;

    double eCurr = eRight;
    double dactCurr = act(Z43*eCurr, lambda) - exact;

    double eNext = 0.0;
    double dactNext = 1.0;

    // double errorRel = std::abs(ePrev - eCurr)/std::abs(ePrev + eCurr + tolerance);
    // double errorAbs = std::abs(ePrev - eCurr);

    while (std::abs(dactCurr - dactPrev) > tolerance && nStep < 100) {

        eNext = eCurr - dactCurr*(eCurr - ePrev)/(dactCurr - dactPrev);
        dactNext = act(Z43*eNext, lambda) - exact;

        ePrev    = eCurr;
        dactPrev = dactCurr;
        eCurr    = eNext;
        dactCurr = dactNext;

        // errorRel = std::abs(ePrev - eCurr)/std::abs(ePrev + eCurr + tolerance);
        // errorAbs = std::abs(ePrev - eCurr);

//        std::cout << "nStep = " << nStep << ", eNext = " << eNext << ", dactNext = " << dactNext << std::endl;

        ++nStep;
    }
    
    eLevel[n][l] = eNext;
    eReady[n][l] = true;
}