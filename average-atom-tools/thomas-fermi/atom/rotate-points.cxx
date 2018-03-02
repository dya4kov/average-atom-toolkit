#include <cmath>
#include <average-atom-tools/thomas-fermi/atom/rotate-points.h>

using namespace AATools::TF;

double RotatePoints::inner(const double& _e, const double& _l) {
    setParam(_e, _l);
    if (!rp1ready) setRP1();
    return rp1;
}

double RotatePoints::outer(const double& _e, const double& _l) {
    setParam(_e, _l);
    if (!rp2ready) setRP2();
    return rp2;
}

Array<RHSRP1::dim> RotatePoints::innerY(const double& _e, const double& _l) {
    setParam(_e, _l);
    if (!rp1ready) setRP1();
    return y1;
}

Array<RHSRP2::dim> RotatePoints::outerY(const double& _e, const double& _l) {
    setParam(_e, _l);
    if (!rp2ready) setRP2();
    return y2;
}

RotatePoints::RotatePoints() :
    V(1.0), T(1.0), Z(1.0),
    e(1.0), l(1.0),
    tolerance(1e-6) 
{
    mu.setZ(1.0);
    mu.setTolerance(1e-10);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhsRP1.set_V(V1);   rhsRP2.set_V(V1);
    rhsRP1.set_T(T1);   rhsRP2.set_T(T1);
    rhsRP1.set_mu(mu1); rhsRP2.set_mu(mu1);

    rp1ready = false; rp2ready = false;
}

RotatePoints::RotatePoints(const RotatePoints& rps) {
    tolerance = rps.tolerance;

    V = rps.V; T = rps.T; Z = rps.Z;

    mu = rps.mu;

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rp1 = rps.rp1; y1 = rps.y1; rp1ready = rps.rp1ready;
    rp2 = rps.rp2; y2 = rps.y2; rp2ready = rps.rp2ready;

    rhsRP1.set_V(V1);   rhsRP2.set_V(V1);
    rhsRP1.set_T(T1);   rhsRP2.set_T(T1);
    rhsRP1.set_mu(mu1); rhsRP2.set_mu(mu1);
}

RotatePoints& RotatePoints::operator=(const RotatePoints& rps) {
    tolerance = rps.tolerance;

    V = rps.V; T = rps.T; Z = rps.Z;

    mu = rps.mu;
    
    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rp1 = rps.rp1; y1 = rps.y1; rp1ready = rps.rp1ready;
    rp2 = rps.rp2; y2 = rps.y2; rp2ready = rps.rp2ready;

    rhsRP1.set_V(V1);   rhsRP2.set_V(V1);
    rhsRP1.set_T(T1);   rhsRP2.set_T(T1);
    rhsRP1.set_mu(mu1); rhsRP2.set_mu(mu1);

    return *this;
}

void RotatePoints::setV(const double& _V) {
    V = _V;

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;

    rhsRP1.set_V(V1);   rhsRP2.set_V(V1);
    rhsRP1.set_mu(mu1); rhsRP2.set_mu(mu1);

    rp1ready = false; rp2ready = false;
}

void RotatePoints::setT(const double& _T) {
    T = _T;

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhsRP1.set_T(T1);   rhsRP2.set_T(T1);
    rhsRP1.set_mu(mu1); rhsRP2.set_mu(mu1);

    rp1ready = false; rp2ready = false;
}

void RotatePoints::setZ(const double& _Z) {
    Z = _Z;

    mu.setZ(Z);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhsRP1.set_V(V1);   rhsRP2.set_V(V1);
    rhsRP1.set_T(T1);   rhsRP2.set_T(T1);
    rhsRP1.set_mu(mu1); rhsRP2.set_mu(mu1);

    rp1ready = false; rp2ready = false;
}

void RotatePoints::setTolerance(const double& t) {
    tolerance = t;
    rp1ready = false; 
    rp2ready = false;
}

void RotatePoints::setParam(const double& _e, const double& _l) {
    if (std::abs(e - _e) > 1e-10 || std::abs(l - _l)) {
        e = _e; l = _l;
        rp1ready = false; 
        rp2ready = false;
    }
}

void RotatePoints::setRP2() {
    y2.fill(0.0);
    rp2ready = true;
    rp2 = 1.0;

    if (e > l) return;

    rhsRP2.reset();
    rhsRP2.set_e(e);
    rhsRP2.set_l(l);

    double t = 1e-10;
    double h = 1e-9;

    double xmin = 0.0;
    double xmax = 1.0;

    Array<RHSRP2::dim> rp2y = y2;

    double error = std::abs(xmax - xmin);

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

    y2  = rhsRP2.yDown();
    rp2 = xmin*xmin;
    return;
}

void RotatePoints::setRP1() {

    if (!rp2ready) setRP2();

    rp1ready = true;

    rhsRP1.reset();
    rhsRP1.set_e(e);
    rhsRP1.set_l(l);

    double t = 1e-10;
    double h = 1e-9;

    double xmin = 0.0;
    double xmax = std::sqrt(rp2);

    Array<RHSRP1::dim> rp1y = y2;

    double error = std::abs(xmax - xmin);

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
    y1  = rhsRP1.yUp();
    rp1 = xmax*xmax;
    return;
}