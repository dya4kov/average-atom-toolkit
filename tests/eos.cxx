#include <average-atom-toolkit/util/unit.h>
#include <average-atom-toolkit/util/quantity.h>
#include <average-atom-toolkit/eos.h>

using namespace aatk;

int main() {
	double Z = 13.0;
	double A = 27.0;

	EoS eos(Z, A);

	Quantity D = 5.0*unit::g/cm3;
	Quantity T = 9.0*unit::eV;

	QuantityArray D(5.0*unit::g/cm3, 15.0*unit::g/cm3, 10);
	QuantityArray T(9.0*unit::eV, 19.0*unit::eV, 10);

	auto data = eos[TF + QE]({P,E})(D, T);

	std::cout << data[P](unit::GPa) << std::endl;

	return 0;
}