#include <cmath>
#include <average-atom-tools/thomas-fermi/atom/action.h>

using namespace AATools::TF;

double Action::operator()(const double& _e, const double& _l) {
    setParam(_e, _l);
    if (!ready) setAction();
    return action;
}

Action::Action() :
    V(1.0), T(1.0), Z(1.0),
    e(1.0), l(1.0),
    tolerance(1e-6) 
{
    mu.setTolerance(1e-10);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    RP.setTolerance(tolerance);

    solver.setTolerance(0.1*tolerance, 0.0);
    ready = false;
}

Action::Action(const Action& a) {
    tolerance = a.tolerance;

    V = a.V; T = a.T;  Z = a.Z; 

    mu = a.mu;

    action = a.action; 
    ready  = a.ready;

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    solver.setTolerance(0.1*tolerance, 0.0);
    
    RP = a.RP;
}

Action& Action::operator=(const Action& a) {
    tolerance = a.tolerance;

    V = a.V; T = a.T;  Z = a.Z; 

    mu = a.mu;

    action = a.action; 
    ready  = a.ready;

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);

    solver.setTolerance(0.1*tolerance, 0.0);
    
    RP = a.RP;

    return *this;
}

void Action::setV(const double& _V) {
    V = _V;
    
    RP.setV(V);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;

    rhs.set_V(V1);
    rhs.set_mu(mu1);
    ready = false;
}

void Action::setT(const double& _T) {
    T = _T;

    RP.setT(T);

    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  T1 = T*std::pow(Z, -4.0/3.0);

    rhs.set_T(T1);
    rhs.set_mu(mu1);
    ready = false;
}

void Action::setZ(const double& _Z) {
    Z = _Z;
    mu.setZ(Z);
    RP.setZ(Z);
    
    double mu1 = mu(V, T)*std::pow(Z, -4.0/3.0);
    double  V1 = V*Z;
    double  T1 = T*std::pow(Z, -4.0/3.0);

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

    double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);
    
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