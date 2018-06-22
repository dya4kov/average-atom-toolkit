#include <cmath>

#include <numeric-toolkit/ODE/types.h>
#include <numeric-toolkit/ODE/solver.h>
#include <numeric-toolkit/ODE/stepper/PD853.h>

#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/ODE/action.h>
#include <average-atom-toolkit/thomas-fermi/atom/action.h>

using namespace aatk::TF;
using numtk::ODE::Array;
using numtk::ODE::Solver;
using numtk::ODE::stepper::PD853;
using ODE::RHSAction;

double Action::operator()(const double& _energy, const double& _lambda) {
    setParam(_energy, _lambda);
    if (!ready) setAction();
    return action;
}

Action::Action() :
    V(1.0), T(1.0), Z(1.0), mu(4.100577730112),
    energy(1.0), lambda(1.0),
    tolerance(1e-6) 
{
    ready = false;
}

Action::Action(const Action& a) {
    tolerance = a.tolerance;
    V = a.V; T = a.T;  Z = a.Z; mu = a.mu;
    energy = a.energy; lambda = a.lambda;
    RP = a.RP;
    action = a.action; 
    ready  = a.ready;
}

Action& Action::operator=(const Action& a) {
    tolerance = a.tolerance;
    V = a.V; T = a.T; Z = a.Z; mu = a.mu;
    energy = a.energy; lambda = a.lambda;
    RP = a.RP;
    action = a.action; 
    ready  = a.ready;
    return *this;
}

void Action::setV(const double& _V) {
    V = _V;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setV(V);
    ready = false;
}

void Action::setT(const double& _T) {
    T = _T;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setT(T);
    ready = false;
}

void Action::setZ(const double& _Z) {
    Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setZ(Z);
    ready = false;
}

void Action::setVTZ(
    const double& _V,
    const double& _T,
    const double& _Z
) {
    V = _V; T = _T; Z = _Z;
    ChemicalPotential M;
    M.setZ(Z);
    M.setTolerance(1e-10);
    mu = M(V, T);
    RP.setVTZ(V, T, Z);
    ready = false;
}

void Action::setTolerance(const double& t) {
    tolerance = t;
    RP.setTolerance(t);
    ready = false;
}

void Action::setParam(const double& _energy, const double& _lambda) {
    if (std::abs(energy - _energy) > 1e-10 || std::abs(lambda - _lambda) > 1e-10) {
        energy = _energy; lambda = _lambda;
        ready = false; 
    }
}

void Action::setAction() {
    double xmax = std::sqrt(RP.outer(energy, lambda));
    double xmin = std::sqrt(RP.inner(energy, lambda));

    ready = true;
    action = 0.0;

    if (xmax <= xmin) return;

    RHSAction rhs;
    Solver<PD853<RHSAction>> solver;
    solver.setTolerance(0.1*tolerance, 0.0);

    Array<RHSAction::dim> ay;
    ay[0] = RP.outerY(energy, lambda)[0]; 
    ay[1] = RP.outerY(energy, lambda)[1];
    ay[2] = 0.0;

    double mu1 = mu*std::pow(Z, -4.0/3.0);
    double V1  = V*Z;
    double T1  = T*std::pow(Z, -4.0/3.0);
    double r0  = std::pow(3.0*V1 / 4.0 / M_PI, 1.0 / 3.0);
    double  l  = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
    double  e  = energy*std::pow(Z, -4.0/3.0);
    
    action = 0.0;

    rhs.set_V(V1);
    rhs.set_T(T1);
    rhs.set_mu(mu1);
    rhs.set_e(e);
    rhs.set_l(l);
    
    solver.setStep(1e-10);
    solver.integrate(rhs, ay, xmax, xmin);

    action = 2.0*r0*std::sqrt(2.0)*ay[2]*std::pow(Z, 1.0/3.0);
    if (action < tolerance) action = 1e+10;

    return;
}