#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom.h>

using namespace aatk::TF;

Atom::Atom() :
    V(1.0), T(1.0), Z(1.0),
    tolerance(1e-6) { }

void Atom::setTolerance(const double& t) {
    tolerance = t;
    e .setTolerance(tolerance);
    RP.setTolerance(tolerance);
    xU.setTolerance(tolerance);
    N .setTolerance(tolerance);
}

void Atom::setV(const double& _V) {
    V = _V;
    e .setV(V);
    RP.setV(V);
    xU.setV(V);
    N .setV(V);
}

void Atom::setT(const double& _T) {
    T = _T;
    e .setT(T);
    RP.setT(T);
    xU.setT(T);
    N .setT(T);
}

void Atom::setZ(const double& _Z) {
    Z = _Z;
    e .setZ(Z);
    RP.setZ(Z);
    xU.setZ(Z);
    N .setZ(Z);
}

void Atom::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    V = _V; T = _T; Z = _Z;
    e .setVTZ(V, T, Z);
    RP.setVTZ(V, T, Z);
    xU.setVTZ(V, T, Z);
    N .setVTZ(V, T, Z);
}

// double Atom::tauEval(const double& e, const double& l) {
// 
//     ::numtk::ODE::Array<RHSTau::dim> y; y.fill(0.0);
//     double rp2 = rPoint2(e, l);
//     double rp1 = rPoint1(e, l);
//     double result = 0.0;
//     double Z = V1/VZ;
// 
//     if (rp2 > rp1) {
//         y[0] = xU(rp2);
//         y[1] = xUDX(rp2);
// 
//         rhsTau.set_e(e);
//         rhsTau.set_l(l);
//         
//         tauSolver.setStep(1e-10);
//         tauSolver.integrate(rhsTau, y, std::sqrt(rp2), std::sqrt(rp1));
// 
//         result = 2.0*std::sqrt(2.0)*r0()*y[2]/Z;
//     }
// 
//     return result;
// }

// double Atom::tau(const int& n, const int& l) {
//     if (needAdjust) adjust();
// 
//     if (tLevel.size() < idxLevel(n + 1) + 1) {
//         tLevel.resize(idxLevel(n + 1) + 1);
//         tReady.resize(idxLevel(n + 1) + 1, false);
//     }
// 
//     double Z = V1/VZ;
// 
//     if (tReady[idxLevel(n) + l]) 
//         return tLevel[idxLevel(n) + l]*r0()*std::pow(Z, -1.0/3.0);
// 
//     double lambda = l + 0.5;
//     double lArg = 0.5*lambda*lambda / r0() / r0() * std::pow(Z, -2.0/3.0);
//     double eArg = e(n,l) * std::pow(Z, -4.0/3.0);
//     tLevel[idxLevel(n) + l] = tauEval(eArg, lArg);
// 
//     return tLevel[idxLevel(n) + l];//*r0()*std::pow(Z, -1.0/3.0);
// }

double Atom::rpInner(const int& n, const int& l) {
    double lambda = l + 0.5;
    double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);
    double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double eArg = std::pow(Z, -4.0/3.0)*e(n, l);
    return RP.inner(eArg, lArg);
}

double Atom::rpOuter(const int& n, const int& l) {
    double lambda = l + 0.5;
    double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);
    double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double eArg = std::pow(Z, -4.0/3.0)*e(n, l);
    return RP.outer(eArg, lArg);
}