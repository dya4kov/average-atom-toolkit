#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>

#include <numeric-tools/ODE/types.h>
#include <average-atom-tools/thomas-fermi/atom.h>

AATools::TF::Atom::Atom() {
    V1 = 1.0; T1 = 1.0;
    VZ = 1.0; TZ = 1.0;
    eps = 1e-6;
    actSolver.setTolerance(0.1*eps, 0.0);
    potSolver.setTolerance(0.0, 0.1*eps);
    tauSolver.setTolerance(0.1*eps, 0.0);
    eStart = {-1e+3, -1e+2, -1e+1, -1.0, 0.0, 1.0, 1e+1, 1e+2, 1e+3, 1e+4};
    needAdjust = true;
}

void AATools::TF::Atom::setTolerance(const double& _eps) {
    eps = _eps;
    actSolver.setTolerance(0.1*eps, 0.0);
    potSolver.setTolerance(0.0, 0.1*eps);
    tauSolver.setTolerance(0.1*eps, 0.0);
}

void AATools::TF::Atom::setV(const double& V) {
    double Z = V1/VZ;
    VZ = V;
    V1 = V*Z;
    needAdjust = true;
}

void AATools::TF::Atom::setT(const double& T) {
    double Z = V1/VZ;
    TZ = T;
    T1 = T*std::pow(Z, -4.0/3.0);
    needAdjust = true;
}

void AATools::TF::Atom::setZ(const double& Z) {
    double Zold = V1/VZ;
    V1 = V1*Z/Zold;
    VZ = V1/Z;
    T1 = T1*std::pow(Z/Zold, -4.0/3.0);
    TZ = T1*std::pow(Z, 4.0/3.0);
    needAdjust = true;
}

void AATools::TF::Atom::setVTZ(const double& V, const double& T, const double& Z) {
    VZ = V;
    TZ = T;
    V1 = V*Z; 
    T1 = T*std::pow(Z, -4.0/3.0);
    needAdjust = true;
}

int AATools::TF::Atom::idxLevel(const int& n) {
    return n*(n - 1)/2;
}

double AATools::TF::Atom::r0() {
    return std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
}

double AATools::TF::Atom::action(const double& e, const double& l) {

    ::numtools::ODE::Array<RHSAction::dim> y; y.fill(0.0);
    double rp2 = rPoint2(e, l);
    double rp1 = rPoint1(e, l);
    double result = 0.0;
    double Z = V1/VZ;

    if (rp2 > rp1) {
        auto rp2y = rhsRP2.yDown();

        y[0] = rp2y[0]; 
        y[1] = rp2y[1];

        rhsAct.set_e(e);
        rhsAct.set_l(l);
        
        actSolver.setStep(1e-10);
        actSolver.integrate(rhsAct, y, std::sqrt(rp2), std::sqrt(rp1));

        result = 2.0*r0()*std::sqrt(2.0)*y[2]*std::pow(Z, 1.0/3.0);
        if (result < eps) result = 1e+10;
    }

    return result;
}

double AATools::TF::Atom::tauEval(const double& e, const double& l) {

    ::numtools::ODE::Array<RHSTau::dim> y; y.fill(0.0);
    double rp2 = rPoint2(e, l);
    double rp1 = rPoint1(e, l);
    double result = 0.0;
    double Z = V1/VZ;

    if (rp2 > rp1) {
        y[0] = xU(rp2);
        y[1] = xUDX(rp2);

        rhsTau.set_e(e);
        rhsTau.set_l(l);
        
        tauSolver.setStep(1e-10);
        tauSolver.integrate(rhsTau, y, std::sqrt(rp2), std::sqrt(rp1));

        result = 2.0*std::sqrt(2.0)*r0()*y[2]/Z;
    }

    return result;
}

double AATools::TF::Atom::rPoint2(const double& e, const double& l) {

    if (e > l) return 1.0;

    ::numtools::ODE::Array<RHSRP2::dim> rp2y; rp2y.fill(0.0);

    rhsRP2.reset();
    rhsRP2.set_e(e);
    rhsRP2.set_l(l);

    double epsNew = 1e-10;
    double h      = 1e-9;

    double xmin = 0.0;
    double xmax = 1.0;

    double error = std::abs(xmax - xmin);

    int nSteps = 0;

    while (error > eps/10 && nSteps < 20) {

        RP2Solver.setStep(h);
        RP2Solver.setTolerance(0.0, epsNew);
        RP2Solver.integrate(rhsRP2, rp2y, xmax, xmin);

        xmin = rhsRP2.xDown();
        xmax = rhsRP2.xUp();

        error = std::abs(xmax - xmin);
        rp2y = rhsRP2.yUp();
        h = error  / 11.0;
        epsNew = h / 21.0;

        ++nSteps;
    }

    return xmin*xmin;
}

double AATools::TF::Atom::rPoint1(const double& e, const double& l) {

    ::numtools::ODE::Array<RHSRP2::dim> rp1y; rp1y.fill(0.0);

    rhsRP1.reset();
    rhsRP1.set_e(e);
    rhsRP1.set_l(l);

    double epsNew = 1e-10;
    double h      = 1e-9;

    double rp2  = rPoint2(e, l);
    double xmin = 0.0;
    double xmax = rhsRP2.xDown();
    rp1y = rhsRP2.yDown();

    double error = std::abs(xmax - xmin);

    int nSteps = 0;

    while (error > eps/10 && nSteps < 20) {

        RP1Solver.setStep(h);
        RP1Solver.setTolerance(0.0, epsNew);
        RP1Solver.integrate(rhsRP1, rp1y, xmax, xmin);

        xmin = rhsRP1.xDown(); 
        xmax = rhsRP1.xUp();

        error = std::abs(xmax - xmin);
        rp1y = rhsRP1.yUp();
        h = error  / 11.0;
        epsNew = h / 21.0;

        ++nSteps;
    }

    return xmax*xmax;
}

double AATools::TF::Atom::e(const int& n, const int& l) {

    if (needAdjust) adjust();
    if (eLevel.size() < idxLevel(n + 1) + 1) {
        eLevel.resize(idxLevel(n + 1) + 1);
        eReady.resize(idxLevel(n + 1) + 1, false);
    }

    if (eReady[idxLevel(n) + l]) 
        return TZ/T1*eLevel[idxLevel(n) + l];
    
    double exact = M_PI*(n - l - 0.5);
    double Z = V1/VZ;
    double lambda = l + 0.5;
    double lArg = 0.5*lambda*lambda / r0() / r0() * std::pow(Z, -2.0/3.0);

    double eLeft;
    double eRight;

    int ie = 0;
    double act = 0.0;
    while (act < exact && ie < eStart.size()) {
        act = action(eStart[ie], lArg);
        ++ie;
    }

    eLeft  = eStart[ie - 2];
    eRight = eStart[ie - 1];

    if (eReady[0]) eLeft = eLevel[0];

    int nSteps = 0;
    double err = std::abs(eLeft - eRight) / std::abs(eLeft + eRight);
    while (err > eps / 10) {
        double eArg = 0.5*(eRight + eLeft);
        double act  = action(eArg, lArg);
        if (act - exact > 0.0) {
            eRight -= 0.5*(eRight - eLeft);
        }
        else {
            eLeft += 0.5*(eRight - eLeft);
        }
        err = std::abs(eLeft - eRight) / std::abs(eLeft + eRight);
        //std::cout << "act = " << act << ", e = " << eArg << std::endl;
        //    std::cout << "rp1 = " << rPoint1(eArg, lArg) << ", rp2 = " << rPoint2(eArg, lArg) << std::endl;
        ++nSteps;
    }

    //if (nSteps == 100 ) {
    //    std::cout << "Warning: energy level convergence failed" << std::endl;
    //}

    //if (std::abs(eMin - eRight) < 10*eps) {
    //    std::cout << "Warning: energy range is small" << std::endl;
    //    std::cout << "Eleft = " << eMin << ", Eright = " << eMax << std::endl;
    //}

    //if (std::abs(eMax - eLeft) < 10*eps) {
    //    std::cout << "Warning: energy range is small" << std::endl;
    //    std::cout << "Eleft = " << eMin << ", Eright = " << eMax << std::endl;
    //}

    eLevel[idxLevel(n) + l] = 0.5*(eLeft + eRight);
    eReady[idxLevel(n) + l] = true;
    return TZ/T1*eLevel[idxLevel(n) + l];
}

double AATools::TF::Atom::rpIn(const int& n, const int& l) {
    if (needAdjust) adjust();

    if (rPointI.size() < idxLevel(n + 1) + 1) {
        rPointI.resize(idxLevel(n + 1) + 1);
        rpIready.resize(idxLevel(n + 1) + 1, false);
    }

    double Z = V1/VZ;

    if (rpIready[idxLevel(n) + l]) 
        return rPointI[idxLevel(n) + l]*r0()* std::pow(Z, -1.0/3.0);

    double lambda = l + 0.5;
    double lArg = 0.5*lambda*lambda / r0() / r0() * std::pow(Z, -2.0/3.0);
    double eArg = e(n,l) * std::pow(Z, -4.0/3.0);
    rPointI[idxLevel(n) + l] = rPoint1(eArg, lArg);

    return rPointI[idxLevel(n) + l];//*r0()*std::pow(Z, -1.0/3.0);
}

double AATools::TF::Atom::rpOut(const int& n, const int& l) {
    if (needAdjust) adjust();

    if (rPointO.size() < idxLevel(n + 1) + 1) {
        rPointO.resize(idxLevel(n + 1) + 1);
        rpOready.resize(idxLevel(n + 1) + 1, false);
    }

    double Z = V1/VZ;

    if (rpOready[idxLevel(n) + l]) 
        return rPointO[idxLevel(n) + l]*r0()*std::pow(Z, -1.0/3.0);

    double lambda = l + 0.5;
    double lArg = 0.5*lambda*lambda / r0() / r0() * std::pow(Z, -2.0/3.0);
    double eArg = e(n,l) * std::pow(Z, -4.0/3.0);
    rPointO[idxLevel(n) + l] = rPoint2(eArg, lArg);

    return rPointO[idxLevel(n) + l];//*r0()*std::pow(Z, -1.0/3.0);
}

double AATools::TF::Atom::tau(const int& n, const int& l) {
    if (needAdjust) adjust();

    if (tLevel.size() < idxLevel(n + 1) + 1) {
        tLevel.resize(idxLevel(n + 1) + 1);
        tReady.resize(idxLevel(n + 1) + 1, false);
    }

    double Z = V1/VZ;

    if (tReady[idxLevel(n) + l]) 
        return tLevel[idxLevel(n) + l]*r0()*std::pow(Z, -1.0/3.0);

    double lambda = l + 0.5;
    double lArg = 0.5*lambda*lambda / r0() / r0() * std::pow(Z, -2.0/3.0);
    double eArg = e(n,l) * std::pow(Z, -4.0/3.0);
    tLevel[idxLevel(n) + l] = tauEval(eArg, lArg);

    return tLevel[idxLevel(n) + l];//*r0()*std::pow(Z, -1.0/3.0);
}

double AATools::TF::Atom::xU(const double& x) {
    if (needAdjust) adjust();
    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    if (xTo < xFrom) potSolver.integrate(rhsPot, phi, xFrom, xTo);
    return phi[0];
}

double AATools::TF::Atom::xUDX(const double& x) {
    if (needAdjust) adjust();
    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    potSolver.integrate(rhsPot, phi, xFrom, xTo);
    if (xTo < xFrom) potSolver.integrate(rhsPot, phi, xFrom, xTo);
    return phi[1];
}

double* AATools::TF::Atom::xU(double* x, size_t n) {
    if (needAdjust) adjust();

    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](size_t i1, size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        potSolver.setStep(eps);
        potSolver.integrate(rhsPot, phi, xFrom, xTo);
        result[i] = phi[0];
        xFrom = xTo;
    }

    return result;
}

double* AATools::TF::Atom::xUDX(double* x, size_t n) {
    if (needAdjust) adjust();

    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](size_t i1, size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        potSolver.setStep(eps);
        potSolver.integrate(rhsPot, phi, xFrom, xTo);
        result[i] = phi[1];
        xFrom = xTo;
    }

    return result;
}

std::vector<double>& AATools::TF::Atom::xU(const std::vector<double>& x) {
    if (needAdjust) adjust();

    size_t n = x.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double>* result = new std::vector<double>(n);

    std::sort(idx.begin(), idx.end(),
       [&x](size_t i1, size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        potSolver.setStep(eps);
        potSolver.integrate(rhsPot, phi, xFrom, xTo);
        (*result)[i] = phi[0];
        xFrom = xTo;
    }

    return *result;
}

std::vector<double>& AATools::TF::Atom::xUDX(const std::vector<double>& x) {
    if (needAdjust) adjust();

    size_t n = x.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double>* result = new std::vector<double>(n);

    std::sort(idx.begin(), idx.end(),
       [&x](size_t i1, size_t i2) {
        return x[i1] > x[i2];
       }
    );

    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        potSolver.setStep(eps);
        potSolver.integrate(rhsPot, phi, xFrom, xTo);
        (*result)[i] = phi[1];
        xFrom = xTo;
    }

    return *result;
}

double AATools::TF::Atom::eDens(const double& x) {
    if (needAdjust) adjust();
    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;
    double xTo   = std::sqrt(x);
    if (xTo < xFrom) potSolver.integrate(rhsPot, phi, xFrom, xTo);
    return TZ*std::sqrt(2.0*TZ)*FDhalf(phi[0]/(x*T1) + mu/T1)/(M_PI*M_PI);
}

double* AATools::TF::Atom::eDens(double* x, size_t n) {
    if (needAdjust) adjust();

    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    double* result = new double[n];

    std::sort(idx.begin(), idx.end(),
       [x](size_t i1, size_t i2) {
        return x[i1] > x[i2]; 
       }
    );

    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);

    double xFrom = 1.0;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        potSolver.setStep(eps);
        potSolver.integrate(rhsPot, phi, xFrom, xTo);
        result[i] = TZ*std::sqrt(2.0*TZ)*FDhalf(phi[0]/(x[i]*T1) + mu/T1)/(M_PI*M_PI);
        xFrom = xTo;
    }

    return result;
}

std::vector<double>& AATools::TF::Atom::eDens(const std::vector<double>& x) {
    if (needAdjust) adjust();

    size_t n = x.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<double>* result = new std::vector<double>(n);

    std::sort(idx.begin(), idx.end(),
       [&x](size_t i1, size_t i2) {
        return x[i1] > x[i2];
       }
    );

    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    double xFrom = 1.0;

    for (auto i : idx) {
        double xTo = std::sqrt(x[i]);
        potSolver.setStep(eps);
        potSolver.integrate(rhsPot, phi, xFrom, xTo);
        (*result)[i] = TZ*std::sqrt(2.0*TZ)*FDhalf(phi[0]/(x[i]*T1) + mu/T1)/(M_PI*M_PI);
        xFrom = xTo;
    }

    return *result;
}

void AATools::TF::Atom::updateMu() {
    ::numtools::ODE::Array<RHSPotential::dim> phi; phi.fill(0.0);
    
    double phi_0 = std::pow(4.0*M_PI/3.0/V1, 1.0/3.0);

    rhsPot.set_V(V1);
    rhsPot.set_T(T1);

    double xFrom = 1.0;
    double xTo   = 0.0;
    double delta = 0.1;
    bool convergenceSuccess = false;

    if (T1 <= 1e-10) T1 = 1e-11;
    double muStart = table(std::log10(V1), std::log10(T1));
    while (!convergenceSuccess) {
        double muLeftStart  = muStart - delta*std::abs(muStart);
        double muRightStart = muStart + delta*std::abs(muStart);
        double muLeft  = muLeftStart;
        double muRight = muRightStart;
        double error   = std::abs(muRight - muLeft)/std::abs(muLeft + muRight + eps);
        while (error > eps) {
            phi.fill(0.0);
            rhsPot.set_mu(0.5*(muLeft + muRight));

            potSolver.integrate(rhsPot, phi, xFrom, xTo);

            if (std::isfinite(phi[0])) {
               if (phi[0] - phi_0 > 0)
                   muRight -= 0.5*(muRight - muLeft);
               else muLeft += 0.5*(muRight - muLeft);
            }
            else {
                muRight -= 0.5*(muRight - muLeft);
            }

            error = std::abs(muRight - muLeft)/std::abs(muRight + muLeft);
        }
        mu = 0.5*(muLeft + muRight);
        convergenceSuccess = true;
        if (std::abs( muLeftStart - muLeft )/std::abs( muLeftStart + muLeft ) < eps*100.0 || 
            std::abs(muRightStart - muRight)/std::abs(muRightStart + muRight) < eps*100.0  ) {
            convergenceSuccess = false;
            delta = delta + 0.1;
            if (delta > 0.45) {
                std::cout << "potential error: too big delta" << std::endl;
                exit(0);
            }
        }
    }
}

void AATools::TF::Atom::adjust() {

    updateMu();

    rhsPot.set_V(V1);
    rhsPot.set_T(T1);
    rhsPot.set_mu(mu);

    rhsAct.set_V(V1); 
    rhsAct.set_T(T1);
    rhsAct.set_mu(mu);

    rhsTau.set_V(V1); 
    rhsTau.set_T(T1);
    rhsTau.set_mu(mu);

    rhsRP1.set_V(V1); 
    rhsRP1.set_T(T1);
    rhsRP1.set_mu(mu);

    rhsRP2.set_V(V1); 
    rhsRP2.set_T(T1);
    rhsRP2.set_mu(mu);

    eLevel.resize(0);
    eReady.resize(0);
    tLevel.resize(0);
    tReady.resize(0);
    rPointI.resize(0);
    rpIready.resize(0);
    rPointO.resize(0);
    rpOready.resize(0);

    needAdjust = false;
}