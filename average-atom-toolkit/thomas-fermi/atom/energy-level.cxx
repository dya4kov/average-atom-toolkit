#include <cmath>
#include <iostream>
#include <thread>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/action.h>

using namespace   aatk::TF;
using ::std::placeholders::_1;
using ::std::placeholders::_2;

EnergyLevel::EnergyLevel() :
    eLevelStart({-1e+3, -1e+2, -1e+1, -1.0, 0.0, 1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5}),
    V(1.0), T(1.0), Z(1.0), threadsLimit(8),
    tolerance(1e-6)
{
    p_runLevel = std::bind(&EnergyLevel::runLevel, this, _1, _2);
    eLevelBuffer.resize(0);
    eLevelReady.resize(0);
}

EnergyLevel::EnergyLevel(const EnergyLevel& e) {
    V = e.V; T = e.T; Z = e.Z;
    action = e.action;
    tolerance = e.tolerance;
    eLevelBuffer = e.eLevelBuffer;
    eLevelReady = e.eLevelReady;
    p_runLevel = std::bind(&EnergyLevel::runLevel, this, _1, _2);
}

EnergyLevel& EnergyLevel::operator=(const EnergyLevel& e) {
    V = e.V; T = e.T; Z = e.Z;
    action = e.action;
    tolerance = e.tolerance;
    eLevelBuffer = e.eLevelBuffer;
    eLevelReady = e.eLevelReady;
    p_runLevel = std::bind(&EnergyLevel::runLevel, this, _1, _2);
    return *this;
}

void EnergyLevel::setV(const double& _V) {
    V = _V;
    action.setV(V);
    eLevelBuffer.resize(0);
    eLevelReady.resize(0);
}

void EnergyLevel::setT(const double& _T) {
    T = _T;
    action.setT(T);
    eLevelBuffer.resize(0);
    eLevelReady.resize(0);
}

void EnergyLevel::setZ(const double& _Z) {
    Z = _Z;
    action.setZ(Z);
    eLevelBuffer.resize(0);
    eLevelReady.resize(0);
}

void EnergyLevel::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    V = _V; T = _T; Z = _Z;
    action.setVTZ(V,T,Z);
    eLevelBuffer.resize(0);
    eLevelReady.resize(0);
}

void EnergyLevel::setThreadsLimit(const std::size_t& Nthreads) {
    threadsLimit = Nthreads; 
}

void EnergyLevel::setTolerance(const double& t) {
    tolerance = t;
    action.setTolerance(t);
    eLevelBuffer.resize(0);
    eLevelReady.resize(0);
}

bool EnergyLevel::bufferOk(const int& n) {
    bool status = false;
    if (eLevelBuffer.size() < iLevel(n + 1) + 1) {
        eLevelBuffer.resize(iLevel(n + 1) + 1);
        eLevelReady .resize(iLevel(n + 1) + 1, false);
    }
    else status = true;
    return status;
}

void EnergyLevel::updateThreads(int& threads, int& current, int& last) {
    for (int thread = current; thread < last; ++thread) {
        if (eLevelReady[thread]) {
            --threads; if (threads == 0) break;
        }
    }
    while (eLevelReady[current] && current < last) ++current;
}

void EnergyLevel::prepareLevelsBelow(const int& nMax) {
    if (bufferOk(nMax)) return;
    int threads = 0;
    int current = 0;
    int last    = 0;
    for (int n = 1; n <= nMax; ++n) {
        for (int l = 0; l < n; ++l) {
            eLevelReady[last] = false;
            std::thread run(p_runLevel, n, l);
            run.detach(); ++threads; ++last;
            while (threads == threadsLimit) {
                updateThreads(threads, current, last);
            }
        }
    }
    bool all_ready = false;
    while (!all_ready) {
        updateThreads(threads, current, last);
        all_ready = (current == last);
    }
    return;
}

double EnergyLevel::operator()(const int& n, const int& l) {
    bufferOk(n);
    if (!eLevelReady[iLevel(n) + l]) {
        runLevel(n, l);
    }
    return std::pow(Z, 4.0/3.0)*eLevelBuffer[iLevel(n) + l];
}

std::vector<int> EnergyLevel::needLevels(const int& n) {
    std::vector<int> list(n);
    list.resize(0);
    for (int l = 0; l < n; ++l)
        if (!eLevelReady[iLevel(n) + l])
            list.push_back(l);
    return list;
}

std::vector<double> EnergyLevel::operator[](const int& n) {
    bufferOk(n);
    auto llist = needLevels(n);

    if (llist.size() > 0) {

        for (auto& l : llist) 
            eLevelReady[iLevel(n) + l] = false;

        int threads = 0;
        for (auto& l : llist) {
            std::thread run(p_runLevel, n, l);
            run.detach(); ++threads;
            while (threads == threadsLimit) {
                for (auto& l : llist) {
                    if (eLevelReady[iLevel(n) + l]) {
                        --threads; if (threads == 0) break;
                    }
                }
            }
        }

        bool all_ready = false;
        while (!all_ready) {
            all_ready = true;
            for (auto& l : llist) 
                all_ready = all_ready && eLevelReady[iLevel(n) + l];
        }
    }

    std::vector<double> result(n);
    for (int l = 0; l < n; ++l)
        result[l] = std::pow(Z, 4.0/3.0)*eLevelBuffer[iLevel(n) + l];

    return result;
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

void EnergyLevel::runLevel(const int& n, const int& l) {
    
    Action act   = action;
    double exact = M_PI*(n - l - 0.5);
    double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);
    double lambda = l + 0.5;
    double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);

    double eLeft;
    double eRight;

    int ie = 0;
    while (act(eLevelStart[ie], lArg) < exact && ie < eLevelStart.size()) ++ie;

    eLeft  = eLevelStart[ie - 1];
    eRight = eLevelStart[ie];

    if (eLevelReady[0]) eLeft = eLevelBuffer[0];

    int nStep = 0;
    double dact = 1.e+5;
    double dactOld = 0.0;
    while (std::abs(dact - dactOld) > 0.1 || std::abs(dact - dactOld) < tolerance) {
        dactOld = dact;
        double eArg = 0.5*(eRight + eLeft);
        dact = act(eArg, lArg) - exact;
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
    double dactPrev = act(ePrev, lArg) - exact;

    double eCurr = eRight;
    double dactCurr = act(eCurr, lArg) - exact;

    double eNext = 0.0;
    double dactNext = 1.0;

    // double errorRel = std::abs(ePrev - eCurr)/std::abs(ePrev + eCurr + tolerance);
    // double errorAbs = std::abs(ePrev - eCurr);

    while (std::abs(dactCurr - dactPrev) > tolerance && nStep < 100) {

        eNext = eCurr - dactCurr*(eCurr - ePrev)/(dactCurr - dactPrev);
        dactNext = act(eNext, lArg) - exact;

        ePrev    = eCurr;
        dactPrev = dactCurr;
        eCurr    = eNext;
        dactCurr = dactNext;

        // errorRel = std::abs(ePrev - eCurr)/std::abs(ePrev + eCurr + tolerance);
        // errorAbs = std::abs(ePrev - eCurr);

//        std::cout << "nStep = " << nStep << ", eNext = " << eNext << ", dactNext = " << dactNext << std::endl;

        ++nStep;
    }

    eLevelBuffer[iLevel(n) + l] = eNext;
    eLevelReady[iLevel(n) + l]  = true;
}