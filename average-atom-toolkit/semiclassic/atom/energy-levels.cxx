#include <cmath>
#include <average-atom-toolkit/semiclassic/atom.h>

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

}
}