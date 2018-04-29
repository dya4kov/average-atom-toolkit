#include <average-atom-toolkit/util/unit.h>
#include <iostream>

using namespace aatk::util::unit;

int main() {

	std::cout << "atomic units: " << std::endl;
	setUnitLength(rBohr);
	setUnitMass(eMass);
	setUnitTime(aTime);

	std::cout << "aVol:    " << aVol.toString() << std::endl;
	std::cout << "eMass:   " << eMass.toString() << std::endl;
	std::cout << "aMass:   " << aMass.toString() << std::endl;
	std::cout << "J:       " << J.toString() << std::endl;
	std::cout << "Ry:      " << Ry.toString() << std::endl;
	std::cout << "hartree: " << hartree.toString() << std::endl;
	std::cout << "eV:      " << eV.toString() << std::endl;
	std::cout << "K:       " << K.toString() << std::endl;
	std::cout << "aTime:   " << aTime.toString() << std::endl;
	std::cout << "hbar:    " << hbar.toString() << std::endl;
	std::cout << "c:       " << c.toString() << std::endl;
	std::cout << "q:       " << qElectron.toString() << std::endl;
	std::cout << "m/s:     " << (m/s).toString() << std::endl;

	std::cout << "SI units: " << std::endl;
	setUnitLength(m);
	setUnitMass(kg);
	setUnitTime(s);

	std::cout << "aVol:    " << aVol.toString() << std::endl;
	std::cout << "eMass:   " << eMass.toString() << std::endl;
	std::cout << "aMass:   " << aMass.toString() << std::endl;
	std::cout << "J:       " << J.toString() << std::endl;
	std::cout << "Ry:      " << Ry.toString() << std::endl;
	std::cout << "hartree: " << hartree.toString() << std::endl;
	std::cout << "eV:      " << eV.toString() << std::endl;
	std::cout << "K:       " << K.toString() << std::endl;
	std::cout << "aTime:   " << aTime.toString() << std::endl;
	std::cout << "hbar:    " << hbar.toString() << std::endl;
	std::cout << "c:       " << c.toString() << std::endl;
	std::cout << "q:       " << qElectron.toString() << std::endl;
	std::cout << "m/s:     " << (m/s).toString() << std::endl;

	std::cout << "SGS units: " << std::endl;
	setUnitLength(cm);
	setUnitMass(g);
	setUnitTime(s);

	std::cout << "aVol:    " << aVol.toString() << std::endl;
	std::cout << "eMass:   " << eMass.toString() << std::endl;
	std::cout << "aMass:   " << aMass.toString() << std::endl;
	std::cout << "J:       " << J.toString() << std::endl;
	std::cout << "Ry:      " << Ry.toString() << std::endl;
	std::cout << "hartree: " << hartree.toString() << std::endl;
	std::cout << "eV:      " << eV.toString() << std::endl;
	std::cout << "K:       " << K.toString() << std::endl;
	std::cout << "aTime:   " << aTime.toString() << std::endl;
	std::cout << "hbar:    " << hbar.toString() << std::endl;
	std::cout << "c:       " << c.toString() << std::endl;
	std::cout << "q:       " << qElectron.toString() << std::endl;
	std::cout << "m/s:     " << (m/s).toString() << std::endl;

	return 0;
}