#include <iostream>
#include <fstream>

#include <average-atom-toolkit/util/quantity-array.h>

using namespace aatk::util;
using namespace aatk::util::unit;

int main() {
	QuantityArray D(5.0*g/cm3, 15.0*g/cm3, 5);
	std::cout << D(g/cm3) << std::endl;

	std::ofstream file("test.txt");
	file << D << std::endl;
	file.close();

	return 0;
}