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
    V(1.0), T(1.0), Z(1.0),
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
}

void EnergyLevel::setT(const double& _T) {
    T = _T;
    action.setT(T);
}

void EnergyLevel::setZ(const double& _Z) {
    Z = _Z;
    action.setZ(Z);
}

void EnergyLevel::setTolerance(const double& t) {
    tolerance = t;
    action.setTolerance(t);
}

void EnergyLevel::checkBufSize(const int& n) {
    if (eLevelBuffer.size() < iLevel(n + 1) + 1) {
        eLevelBuffer.resize(iLevel(n + 1) + 1);
        eLevelReady .resize(iLevel(n + 1) + 1, false);
    }
}

double EnergyLevel::operator()(const int& n, const int& l) {
    checkBufSize(n);
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
    checkBufSize(n);
    auto llist = needLevels(n);

    if (llist.size() > 0) {
        std::vector<std::thread> multi;
        for (auto& l : llist)
            multi.push_back(std::thread(p_runLevel, n, l));    
        for (auto& l : llist)
            multi[l].join();
    }

    std::vector<double> result(n);
    for (int l = 0; l < n; ++l)
        result[l] = std::pow(Z, 4.0/3.0)*eLevelBuffer[iLevel(n) + l];

    return result;
}

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
    double err = std::abs(eLeft - eRight) / std::abs(eLeft + eRight);
    while (err > 0.1*tolerance && nStep < 100) {
        double eArg = 0.5*(eRight + eLeft);
        if (act(eArg, lArg) - exact > 0.0) {
            eRight -= 0.5*(eRight - eLeft);
        }
        else {
            eLeft += 0.5*(eRight - eLeft);
        }
        err = std::abs(eLeft - eRight) / std::abs(eLeft + eRight);
        ++nStep;
    }

    // if (nStep == 100 ) {
    //     std::cout << "Warning: energy level convergence failed" << std::endl;
    // }

    //if (std::abs(eMin - eRight) < 10*tolerance) {
    //    std::cout << "Warning: energy range is small" << std::endl;
    //    std::cout << "Eleft = " << eMin << ", Eright = " << eMax << std::endl;
    //}

    //if (std::abs(eMax - eLeft) < 10*tolerance) {
    //    std::cout << "Warning: energy range is small" << std::endl;
    //    std::cout << "Eleft = " << eMin << ", Eright = " << eMax << std::endl;
    //}

    eLevelBuffer[iLevel(n) + l] = 0.5*(eLeft + eRight);
    eLevelReady[iLevel(n) + l]  = true;
}