#include <iostream>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

using namespace numtools::specfunc;

int main() {

	FermiDirac<FD::DMHalf>    FD_DMHalf;
	FermiDirac<FD::MHalf>     FD_MHalf;
	FermiDirac<FD::Half>      FD_Half;
	FermiDirac<FD::ThreeHalf> FD_3Half;

	std::cout << "I'_{-1/2}(1)  = " << FD_DMHalf(1.0) << std::endl;
	std::cout << "I_{-1/2}(1)   = " << FD_MHalf(1.0)  << std::endl;
	std::cout << "I_{1/2}(1)    = " << FD_Half(1.0)   << std::endl;
	std::cout << "I_{3/2}(1)    = " << FD_3Half(1.0)  << std::endl;

	std::cout << "I'_{-1/2}(10) = " << FD_DMHalf(10.0) << std::endl;
	std::cout << "I_{-1/2}(10)  = " << FD_MHalf(10.0)  << std::endl;
	std::cout << "I_{1/2}(10)   = " << FD_Half(10.0)   << std::endl;
	std::cout << "I_{3/2}(10)   = " << FD_3Half(10.0)  << std::endl;

	return 0;
}