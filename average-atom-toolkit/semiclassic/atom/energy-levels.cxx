#include <cmath>
#include <average-atom-toolkit/semiclassic/atom.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>

extern "C" {
#include <gsl/gsl_integration.h>
}

namespace aatk {
namespace semiclassic {

double Atom::energyLevel(int n, int l) {
	if (!eLevelReady[n][l]) evaluateEnergyLevel(n, l);
	return eLevel[n][l];
}

int Atom::evaluateEnergyLevel(int n, int l) {
    double exact = M_PI*(n - l - 0.5);
    double lambda = l + 0.5;
    int ie = 0;
    while (action(eLevelStart[ie], lambda) < exact && ie < eLevelStart.size()) ++ie;

    double eLeft  = eLevelStart[ie - 1];
    double eRight = eLevelStart[ie];

    int nStep = 0;
    double dact = 1.e+5;
    double dactOld = 0.0;
    while (std::abs(dact - dactOld) > 0.1 || std::abs(dact - dactOld) < tolerance) {
        dactOld = dact;
        double eArg = 0.5*(eRight + eLeft);
        dact = action(eArg, lambda) - exact;
        if (dact > 0.0) {
            eRight -= 0.5*(eRight - eLeft);
        }
        else {
            eLeft += 0.5*(eRight - eLeft);
        }
        ++nStep;
    }

    double ePrev = eLeft;
    double dactPrev = action(ePrev, lambda) - exact;

    double eCurr = eRight;
    double dactCurr = action(eCurr, lambda) - exact;

    double eNext = 0.0;
    double dactNext = 1.0;

    while (std::abs(dactCurr - dactPrev) > tolerance && nStep < 100) {

        eNext = eCurr - dactCurr*(eCurr - ePrev)/(dactCurr - dactPrev);
        dactNext = action(eNext, lambda) - exact;

        ePrev    = eCurr;
        dactPrev = dactCurr;
        eCurr    = eNext;
        dactCurr = dactNext;

        ++nStep;
    }
    
    eLevel[n][l] = eNext;
    eLevelReady[n][l] = true;
    return 0;
}

double Atom::energyDensityContinuous(double x){
    numtk::specfunc::FermiDirac<numtk::specfunc::FD::ThreeHalf>        FD_ThreeHalf;
    numtk::specfunc::FermiDiracInc<numtk::specfunc::FDI::ThreeHalfInc> FD_ThreeHalf_Inc;
    const double E0 = energyLevel(nmax,nmax - 1);
    double V_r ;
    potential(&x,&V_r,1);
    const double mu = chemPot;
    const double T = temperature;
    const double factor = T * pow(2 * T,3.0 / 2.0) / (2 * pow(M_PI, 2) );
    const double y0 = (V_r + E0)/T;
    double result;

    if (y0 < 0){
        result =  FD_ThreeHalf((V_r + mu)/T);
    }
    else{
        result = FD_ThreeHalf((V_r + mu)/T) - FD_ThreeHalf_Inc((V_r + mu)/T,y0); // fermi_dirac_1/2 + fermi_dirac_1/2_incomplete
    }

    return result * factor;
}

double Atom::energyDensityContinuousFunc(double x, void * atomClass){
    Atom * temp_cell = (Atom *)atomClass;

    return  x * x * temp_cell->energyDensityContinuous(x);
}

double Atom::energyContinuous() {
    double result, error;
    Atom * tempAtom = (Atom *)this;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function Func;
    Func.function = energyDensityContinuousFunc;
    Func.params = tempAtom;

    gsl_integration_qags (&Func,1e-6,1,tolerance ,tolerance,1000,w,&result, &error); // epsabs = 1e-6
    gsl_integration_workspace_free (w);

    return 4 * M_PI * result * pow(r0,3);
}

double Atom::energyFull(){
    double result = 0.0;

    for (int n = 1; n < nmax; n++){
        for (int l = 0; l < n; l++ ){
            result += energyLevel(n,l) * electronStates(n, l);
        }
    }

    if (useContinuous){
        result += energyContinuous();
    }

    return result;
}


}
}