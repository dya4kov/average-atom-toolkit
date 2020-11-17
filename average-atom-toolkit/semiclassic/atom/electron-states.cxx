#include <cmath>
#include <average-atom-toolkit/semiclassic/atom.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>
#include <iomanip>

extern "C" {
#include <gsl/gsl_integration.h>
}

namespace aatk {
namespace semiclassic {

double Atom::electronStates(int n, int l) {
	if (!chemPotReady) evaluateChemicalPotential();
	double enl = energyLevel(n, l);
    double exponent = (enl - chemPot)/temperature;
    double Nnl = 0.0;
    if (exponent > 50.0) Nnl = 0.0;
    else if (exponent < -50.0) Nnl = 2.0*(2.0 * l + 1.0); 
    else Nnl = 2.0*(2.0 * l + 1.0) / (1.0 + std::exp(exponent));
    return Nnl;
}
double Atom::electronStates(int n) {
	double Nn = 0.0;
	for (int l = 0; l < n; ++l) Nn += electronStates(n, l);
	return Nn;
}
double Atom::electronStates() {
	double N = 0.0;
	for (int n = 1; n <= nmax; ++n) {
        N += electronStates(n);
    }
	return N;
}

double Atom::electronStates(double CP) {
	double N = 0.0;
	for (int n = 1; n <= nmax; ++n) {
        for (int l = 0; l < n; ++l) {
            double enl = energyLevel(n, l);
            double exponent = (enl - CP)/temperature;
            double Nnl;
            if (exponent > 50.0) Nnl = 0.0;
            else if (exponent < -50.0) Nnl = 2.0*(2.0 * l + 1.0); 
            else Nnl = 2.0*(2.0 * l + 1.0) / (1.0 + std::exp(exponent));
            N += Nnl;
        }
    }
	return N;
}

double Atom::electronDensityContinious(double x) {
    numtk::specfunc::FermiDirac<numtk::specfunc::FD::Half>        FD_Half;
    numtk::specfunc::FermiDiracInc<numtk::specfunc::FDI::HalfInc> FD_Half_Inc;
    const double E_hartre = 27.21;

    const double E0 = energyLevel(nmax,nmax - 1);
    double V_r ;
    potential(&x,&V_r,1);
    const double mu = chemPot;
    const double T = temperature;
    const double factor = pow(2 * T,3.0 / 2.0) / (2 * pow(M_PI, 2) );
    const double y0 = (V_r + E0)/T;
    double result = 0.0;

    if (y0 < 0){
        result =  FD_Half((V_r + mu)/T);
    }
    else{
        result = FD_Half((V_r + mu)/T) - FD_Half_Inc((V_r + mu)/T,y0); // fermi_dirac_1/2 + fermi_dirac_1/2_incomplete
    }

    return result * factor;
}

double Atom::electronStatesContinuousFunc (double x, void * atomClass){
    Atom * temp_cell = (Atom *)atomClass;
    return  x * x * temp_cell->electronDensityContinious(x);
}

double Atom::electronStatesContinuous(double CP){
    double result, error;
    double oldChemPot = chemPot;
    chemPot = CP;
    Atom * tempAtom = (Atom *)this;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function Func;
    Func.function = electronStatesContinuousFunc;
    Func.params = tempAtom;

    gsl_integration_qags (&Func,1e-6,1,0.0 ,tolerance,1000,w,&result, &error);
    gsl_integration_workspace_free (w);

    chemPot = oldChemPot;

    return 4 * M_PI * result * pow(r0,3);
}


void Atom::evaluateChemicalPotential() {
	double CPprev = chemPot - 0.25*std::abs(chemPot);
    double Nprev = electronStates(CPprev) - Zcharge;
    double CPcurr = chemPot + 0.25*std::abs(chemPot);
    double Ncurr = electronStates(CPcurr) - Zcharge;

    if(useContinuous){
        Nprev += electronStatesContinuous(CPprev);
        Ncurr += electronStatesContinuous(CPcurr);
    }

    double CPnext = 0.0;
    double Nnext = 1.0;

    int nStep = 0;

    while (std::abs(Ncurr - Nprev) > tolerance && nStep < 100) {

        CPnext = CPcurr - Ncurr*(CPcurr - CPprev)/(Ncurr - Nprev);
        Nnext = electronStates(CPnext) - Zcharge;

        if(useContinuous){
            Nnext += electronStatesContinuous(CPnext);
        }

        CPprev = CPcurr;
        Nprev  = Ncurr;
        CPcurr = CPnext;
        Ncurr  = Nnext;

        ++nStep;
    }

    chemPot = CPnext;
    chemPotReady = true;
}

}
}