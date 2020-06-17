#include <cmath>
#include <average-atom-toolkit/semiclassic/atom.h>

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

void Atom::evaluateChemicalPotential() {
	double CPprev = chemPot - 0.25*std::abs(chemPot);
    double Nprev = electronStates(CPprev) - Zcharge;

    double CPcurr = chemPot + 0.25*std::abs(chemPot);
    double Ncurr = electronStates(CPcurr) - Zcharge;

    double CPnext = 0.0;
    double Nnext = 1.0;

    int nStep = 0;

    while (std::abs(Ncurr - Nprev) > tolerance && nStep < 100) {

        CPnext = CPcurr - Ncurr*(CPcurr - CPprev)/(Ncurr - Nprev);
        Nnext = electronStates(CPnext) - Zcharge;

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