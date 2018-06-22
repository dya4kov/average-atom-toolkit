#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>
#include <average-atom-toolkit/thomas-fermi/atom/ODE/rotate-points.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

using namespace aatk::TF;

using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;
using ODE::RHSRP1;
using ODE::RHSRP2;

double RotatePoints::inner(const double& _energy, const double& _lambda) {
    setParam(_energy, _lambda);
    if (!rpIready) setRPinner();
    return rpInner;
}

double RotatePoints::outer(const double& _energy, const double& _lambda) {
    setParam(_energy, _lambda);
    if (!rpOready) setRPouter();
    return rpOuter;
}

double* RotatePoints::innerY(const double& _energy, const double& _lambda) {
    setParam(_energy, _lambda);
    if (!rpIready) setRPinner();
    return yInner;
}

double* RotatePoints::outerY(const double& _energy, const double& _lambda) {
    setParam(_energy, _lambda);
    if (!rpOready) setRPouter();
    return yOuter;
}

RotatePoints::RotatePoints() :
    V(1.0), T(1.0), Z(1.0), mu(4.100577730112),
    energy(1.0), lambda(1.0),
    tolerance(1e-6) 
{
    rpIready = false; rpOready = false;
}

RotatePoints::RotatePoints(const RotatePoints& rps) {
    tolerance = rps.tolerance;

    V = rps.V; T = rps.T; Z = rps.Z; mu = rps.mu;
    energy = rps.energy; lambda = rps.lambda;

    rpInner = rps.rpInner; rpIready = rps.rpIready;
    rpOuter = rps.rpOuter; rpOready = rps.rpOready;

    yOuter[0] = rps.yOuter[0];
    yOuter[1] = rps.yOuter[1];
}

RotatePoints& RotatePoints::operator=(const RotatePoints& rps) {
    tolerance = rps.tolerance;

    V = rps.V; T = rps.T; Z = rps.Z; mu = rps.mu;
    energy = rps.energy; lambda = rps.lambda;

    rpInner = rps.rpInner; rpIready = rps.rpIready;
    rpOuter = rps.rpOuter; rpOready = rps.rpOready;

    yOuter[0] = rps.yOuter[0];
    yOuter[1] = rps.yOuter[1];

    return *this;
}

void RotatePoints::setV(const double& _V) {
    V = _V;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    rpIready = false; rpOready = false;
}

void RotatePoints::setT(const double& _T) {
    T = _T;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    rpIready = false; rpOready = false;
}

void RotatePoints::setZ(const double& _Z) {
    Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    rpIready = false; rpOready = false;
}

void RotatePoints::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    V = _V; T = _T; Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    rpIready = false; rpOready = false;
}

void RotatePoints::setTolerance(const double& t) {
    tolerance = t;
    rpIready = false; 
    rpOready = false;
}

void RotatePoints::setParam(const double& _energy, const double& _lambda) {
    if (std::abs(energy - _energy) > 1e-10 || std::abs(lambda - _lambda) > 1e-10) {
        energy = _energy; lambda = _lambda;
        rpIready = false; 
        rpOready = false;
    }
}

void RotatePoints::setRPouter() {
    yOuter[0] = 0.0;
    yOuter[1] = 0.0;
    rpOready = true;
    rpOuter  = 1.0;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    double  r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double   l = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double   e = energy*std::pow(Z, -4.0/3.0);

    if (e > l) return;

    RHSRP2 rhsRP2;
    rhsRP2.reset();
    rhsRP2.set_V(V1);
    rhsRP2.set_T(T1);
    rhsRP2.set_mu(mu1);
    rhsRP2.set_e(e);
    rhsRP2.set_l(l);

    double t = 1e-10;
    double h = 1e-9;

    double xmin = 0.0;
    double xmax = 1.0;

    Array<RHSRP2::dim> rp2y; rp2y.fill(0.0);

    double error = std::abs(xmax - xmin);
    Solver<PD853<RHSRP2>> solverRP2;
    int nStep = 0;

    while (error > 0.1*tolerance && nStep < 20) {

        solverRP2.setStep(h);
        solverRP2.setTolerance(0.0, t);
        solverRP2.integrate(rhsRP2, rp2y, xmax, xmin);

        xmin = rhsRP2.xDown();
        xmax = rhsRP2.xUp();

        error = std::abs(xmax - xmin);
        rp2y  = rhsRP2.yUp();
        h     = error  / 11.0;
        t     = h / 21.0;

        ++nStep;
    }

    yOuter[0] = rhsRP2.yDown()[0];
    yOuter[1] = rhsRP2.yDown()[1];
    rpOuter = xmin*xmin;
    return;
}

void RotatePoints::setRPinner() {

    if (!rpOready) setRPouter();

    rpIready = true;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    double  r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double   l = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double   e = energy*std::pow(Z, -4.0/3.0);

    RHSRP1 rhsRP1;
    rhsRP1.reset();
    rhsRP1.set_V(V1);
    rhsRP1.set_T(T1);
    rhsRP1.set_mu(mu1);
    rhsRP1.set_e(e);
    rhsRP1.set_l(l);

    double t = 1e-10;
    double h = 1e-9;

    double xmin = 0.0;
    double xmax = std::sqrt(rpOuter);

    Array<RHSRP1::dim> rp1y;
    rp1y[0] = yOuter[0];
    rp1y[1] = yOuter[1];

    double error = std::abs(xmax - xmin);
    Solver<PD853<RHSRP1>> solverRP1;
    int nStep = 0;

    while (error > 0.1*tolerance && nStep < 20) {

        solverRP1.setStep(h);
        solverRP1.setTolerance(0.0, t);
        solverRP1.integrate(rhsRP1, rp1y, xmax, xmin);

        xmin = rhsRP1.xDown();
        xmax = rhsRP1.xUp();

        error = std::abs(xmax - xmin);
        rp1y  = rhsRP1.yUp();
        h     = error / 11.0;
        t     = h / 21.0;

        ++nStep;
    }

    yInner[0] = rhsRP1.yUp()[0];
    yInner[1] = rhsRP1.yUp()[1];
    rpInner = xmax*xmax;
    return;
}