#include <cmath>
#include <average-atom-tools/thomas-fermi/atom/action.h>

using namespace AATools::TF;

double Action::operator()(const double& _e, const double& _l) {
    setParam(_e, _l);
    if (!ready) setAction();
    return action;
}

Action::Action() :
    V1(1.0), T1(1.0), e(1.0),
    VZ(1.0), TZ(1.0), l(1.0),
    tolerance(1e-6) 
{
    potential.setTolerance(1e-10);
    potential.setV(VZ);
    potential.setT(TZ);
    potential.setZ(1.0);
    muZ = potential.mu();
    mu1 = muZ;

    RP.setV(V1);  rhs.set_V(V1);
    RP.setT(T1);  rhs.set_T(T1);
    RP.setZ(1.0); rhs.set_mu(mu1);
    RP.setTolerance(tolerance);

    solver.setTolerance(0.1*tolerance, 0.0);
    ready = false;
}

Action::Action(const Action& a) {
    tolerance = a.tolerance;

    V1  = a.V1;  VZ = a.VZ; 
    T1  = a.T1;  TZ = a.TZ; 
    mu1 = a.mu1; muZ = a.muZ;

    action = a.action; 
    ready  = a.ready;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    solver.setTolerance(0.1*tolerance, 0.0);

    potential = a.potential;
    RP        = a.RP;
}

Action& Action::operator=(const Action& a) {
    tolerance = a.tolerance;

    V1  = a.V1;  VZ = a.VZ; 
    T1  = a.T1;  TZ = a.TZ; 
    mu1 = a.mu1; muZ = a.muZ;

    action = a.action; 
    ready  = a.ready;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    solver.setTolerance(0.1*tolerance, 0.0);

    potential = a.potential;
    RP        = a.RP;

    return *this;
}

void Action::setV(const double& V) {
    double Z = V1/VZ;
    VZ = V;
    V1 = V*Z;
    potential.setV(V);
    muZ = potential.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);

    RP.setV(V); 

    rhs.set_V(V1);
    rhs.set_mu(mu1);
    ready = false;
}

void Action::setT(const double& T) {
    double Z = V1/VZ;
    TZ = T;
    T1 = T*std::pow(Z, -4.0/3.0);
    potential.setT(T);
    muZ = potential.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);

    RP.setT(T); 

    rhs.set_T(T1);
    rhs.set_mu(mu1);
    ready = false;
}

void Action::setZ(const double& Z) {
    double Zold = V1/VZ;
    V1 = V1*Z/Zold;
    VZ = V1/Z;
    T1 = T1*std::pow(Z/Zold, -4.0/3.0);
    TZ = T1*std::pow(Z, 4.0/3.0);
    potential.setZ(Z);
    muZ = potential.mu();
    mu1 = muZ*std::pow(Z, -4.0/3.0);

    RP.setZ(Z);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    ready = false;
}

void Action::setTolerance(const double& t) {
    tolerance = t;
    ready = false;
}

void Action::setParam(const double& _e, const double& _l) {
    if (std::abs(e - _e) > 1e-10 || std::abs(l - _l)) {
        e = _e; l = _l;
        ready = false; 
    }
}

void Action::setAction() {
    double xmax = std::sqrt(RP.outer(e, l));
    double xmin = std::sqrt(RP.inner(e, l));

    Array<RHSAction::dim> ay;
    ay[0] = RP.outerY(e, l)[0]; 
    ay[1] = RP.outerY(e, l)[1];
    ay[2] = 0.0;

    double Z = V1/VZ;
    double r0 = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    
    action = 0.0;

    if (xmax > xmin) {

        rhs.set_e(e);
        rhs.set_l(l);
        
        solver.setStep(1e-10);
        solver.integrate(rhs, ay, xmax, xmin);

        action = 2.0*r0*std::sqrt(2.0)*ay[2]*std::pow(Z, 1.0/3.0);
        if (action < tolerance) action = 1e+10;
    }

    return;
}