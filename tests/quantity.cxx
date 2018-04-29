#include <average-atom-toolkit/util/quantity.h>
#include <iostream>

using namespace aatk::util;
using namespace aatk::util::unit;

int main() {
	Quantity D = aMass/aVol;
	std::cout << D(g/cm3) << std::endl;

	return 0;
}