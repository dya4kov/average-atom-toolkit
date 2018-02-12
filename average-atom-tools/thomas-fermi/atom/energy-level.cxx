#include <cmath>
#include <iostream>
#include <thread>

#include <average-atom-tools/thomas-fermi/atom/energy-level.h>
#include <average-atom-tools/thomas-fermi/atom/action.h>

using namespace   AATools::TF;
using ::std::placeholders::_1;
using ::std::placeholders::_2;

EnergyLevel::EnergyLevel() :
    eLevelStart({-1e+3, -1e+2, -1e+1, -1.0, 0.0, 1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5}),
    V1(1.0), T1(1.0),
    VZ(1.0), TZ(1.0),
    tolerance(1e-6)
{
    p_runLevel = std::bind(&EnergyLevel::runLevel, this, _1, _2);
    eLevelBuffer.resize(0);
    eLevelReady.resize(0);
}

EnergyLevel::EnergyLevel(const EnergyLevel& e) {
    V1 = e.V1; VZ = e.VZ;
    T1 = e.T1; TZ = e.TZ;
    tolerance = e.tolerance;
    eLevelBuffer = e.eLevelBuffer;
    eLevelReady = e.eLevelReady;
    p_runLevel = std::bind(&EnergyLevel::runLevel, this, _1, _2);
}

EnergyLevel& EnergyLevel::operator=(const EnergyLevel& e) {
    V1 = e.V1; VZ = e.VZ;
    T1 = e.T1; TZ = e.TZ;
    tolerance = e.tolerance;
    eLevelBuffer = e.eLevelBuffer;
    eLevelReady = e.eLevelReady;
    return *this;
}

void EnergyLevel::setV(const double& V) {
    double Z = V1/VZ;
    VZ = V;
    V1 = V*Z;
}

void EnergyLevel::setT(const double& T) {
    double Z = V1/VZ;
    TZ = T;
    T1 = T*std::pow(Z, -4.0/3.0);
}

void EnergyLevel::setZ(const double& Z) {
    double Zold = V1/VZ;
    V1 = V1*Z/Zold;
    VZ = V1/Z;
    T1 = T1*std::pow(Z/Zold, -4.0/3.0);
    TZ = T1*std::pow(Z, 4.0/3.0);
}

void EnergyLevel::setTolerance(const double& t) {
    tolerance = t;
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
    return TZ/T1*eLevelBuffer[iLevel(n) + l];
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
        result[l] = TZ/T1*eLevelBuffer[iLevel(n) + l];

    return result;
}

void EnergyLevel::runLevel(const int& n, const int& l) {

    double Z = V1/VZ;

    Action action;
    action.setTolerance(tolerance);
    action.setV(VZ);
    action.setT(TZ);
    action.setZ(Z);
    
    double exact = M_PI*(n - l - 0.5);
    double r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double lambda = l + 0.5;
    double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);

    double eLeft;
    double eRight;

    int ie = 0;
    while (action(eLevelStart[ie], lArg) < exact && ie < eLevelStart.size()) ++ie;

    eLeft  = eLevelStart[ie - 1];
    eRight = eLevelStart[ie];

    if (eLevelReady[0]) eLeft = eLevelBuffer[0];

    int nStep = 0;
    double err = std::abs(eLeft - eRight) / std::abs(eLeft + eRight);
    while (err > 0.1*tolerance && nStep < 100) {
        double eArg = 0.5*(eRight + eLeft);
        double act  = action(eArg, lArg);
        if (act - exact > 0.0) {
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