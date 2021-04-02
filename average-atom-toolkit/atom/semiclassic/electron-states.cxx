#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

#include <average-atom-toolkit/atom/semiclassic.h>

namespace aatk {
namespace atom {

double SemiclassicAtom::electronStatesDiscrete(int n, int l) {
	if (!chemPotReady) evaluate_chemical_potential();
	double enl = energyLevel(n, l);
    double exponent = (enl - M)/T;
    double Nnl = 0.0;

    if (exponent > 50.0) Nnl = 0.0;
    else if (exponent < -50.0) Nnl = 2.0*(2.0 * l + 1.0);
    else Nnl = 2.0*(2.0 * l + 1.0) / (1.0 + std::exp(exponent));

    if (useContinuous) {
        Nnl *= (enl < boundaryEnergy) ? 1.0 : 0.0;
    }

    return Nnl;
}
double SemiclassicAtom::electronStatesDiscrete(int n) {
	double Nn = 0.0;
	for (int l = 0; l < n; ++l) Nn += electronStatesDiscrete(n, l);
	return Nn;
}
double SemiclassicAtom::electronStatesDiscrete() {
	double N = 0.0;
	for (int n = 1; n <= nmax; ++n) {
        N += electronStatesDiscrete(n);
    }
	return N;
}
double SemiclassicAtom::electronStatesDiscrete(double CP) {
	double N = 0.0;
	for (int n = 1; n <= nmax; ++n) {
        for (int l = 0; l < n; ++l) {
            double enl = energyLevel(n, l);
            double exponent = (enl - CP)/T;
            double Nnl;

            if (exponent > 50.0) Nnl = 0.0;
            else if (exponent < -50.0) Nnl = 2.0 * (2.0 * l + 1.0);
            else Nnl = 2.0 * (2.0 * l + 1.0) / (1.0 + std::exp(exponent));

            if (useContinuous) {
            Nnl *= (enl < boundaryEnergy) ? 1.0 : 0.0;
            }
                N += Nnl;
        }
    }
	return N;
}

double electronStatesContinuousFunc (double x, void * params){
    SemiclassicAtom * atom = (SemiclassicAtom  *)params;
    return  x * x * atom->electronDensityContinuous(x);
}
double SemiclassicAtom::electronStatesContinuous(double CP){
    double result, error;
    double old_M = M;
    M = CP;
    SemiclassicAtom  * atom = (SemiclassicAtom  *)this;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function Func;

    Func.function = &electronStatesContinuousFunc;
    Func.params = atom;

    gsl_integration_qags (&Func,1e-6,1,tolerance ,tolerance,1000,w,&result, &error);
    gsl_integration_workspace_free (w);
    M = old_M;

    return 4 * M_PI * result * std::pow(r0,3);
}
double SemiclassicAtom::electronStatesContinuous(){
    return electronStatesContinuous(M);
}

struct SCstatesParams {
    SemiclassicAtom *atom;
};

double SCstatesFunc(double x, void *params) {
    auto& p = *((SCstatesParams *) params);
    auto& atom = *(p.atom);
    double result = 0.0;

    if (atom.boundaryEnergyValue() != 0.0 ){
        result = atom.electronStatesDiscrete(x) + atom.electronStatesContinuous(x)
               - atom.Znucleus();
    }
    else{
        result = atom.electronStatesDiscrete(x) - atom.Znucleus();
    }

    return result;
}

void SemiclassicAtom::evaluate_chemical_potential() {

    SCstatesParams params;
    params.atom = this;

    double MMin = energyLevel(1, 0) - 40*T; // ?
    double MMax = 40*T;

    gsl_function F;
    F.function = &SCstatesFunc;
    F.params = &params;

    int status;
    int iter, max_iter = 100;

    // locate root for chemical potential
    gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, MMin, MMax);
    iter = 0;
    status = GSL_CONTINUE;
    while (status == GSL_CONTINUE && iter < max_iter) {
        status = gsl_root_fsolver_iterate (solver);
        MMin   = gsl_root_fsolver_x_lower (solver);
        MMax   = gsl_root_fsolver_x_upper (solver);
        status = gsl_root_test_interval(MMin, MMax, 0.0, tolerance);
    }
    gsl_root_fsolver_free(solver);

    M = 0.5*(MMin + MMax);
    chemPotReady = true;
}

}
}