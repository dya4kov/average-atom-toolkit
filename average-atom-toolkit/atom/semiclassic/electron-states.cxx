#include <cmath>
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
            else if (exponent < -50.0) Nnl = 2.0*(2.0 * l + 1.0); 
            else Nnl = 2.0*(2.0 * l + 1.0) / (1.0 + std::exp(exponent));
            N += Nnl;
        }
    }
	return N;
}

void SemiclassicAtom::evaluate_chemical_potential() {
    double Mprev = M - 0.25*std::abs(M);
    double Nprev = electronStatesDiscrete(Mprev);
    double Mcurr = M + 0.25*std::abs(M);
    double Ncurr = electronStatesDiscrete(Mcurr);

    // if(useContinuous){
    //     Nprev += electronStatesContinuous(CPprev);
    //     Ncurr += electronStatesContinuous(CPcurr);
    // }

    double Mnext = 0.0;
    double Nnext = 1.0;

    int nStep = 0;

    while (std::abs(Ncurr - Nprev) > tolerance && nStep < 100) {

        Mnext = Mcurr - Ncurr*(Mcurr - Mprev)/(Ncurr - Nprev);
        Nnext = electronStatesDiscrete(Mnext) - Z;

        // if(useContinuous){
        //     Nnext += electronStatesContinuous(CPnext);
        // }

        Mprev = Mcurr;
        Nprev = Ncurr;
        Mcurr = Mnext;
        Ncurr = Nnext;

        ++nStep;
    }

    M = Mnext;
    chemPotReady = true;
}

}
}