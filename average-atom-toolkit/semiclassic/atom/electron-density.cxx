#include <cmath>
#include <average-atom-toolkit/semiclassic/atom.h>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>

namespace aatk {
namespace semiclassic {

std::vector<double> Atom::electronDensity(const std::vector<double>& x) {
	std::vector<double> result(x.size(), 0.0);
	electronDensity(x.data(), result.data(), x.size());
	return result;
}
double Atom::electronDensity(double x) {
	double result;
	electronDensity(&x, &result, 1);
	return result;
}
void Atom::electronDensity(const double* x, double* y, std::size_t n) {
	auto& density = *densityInterpolation;
	for (std::size_t i = 0; i < n; ++i) {
		y[i] = density(std::sqrt(x[i]));
	}
}

double Atom::electronDensityContinuous(double x) {
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

void Atom::electronDensityContinuous(const double* x, double* y, std::size_t n) {
    // to do
}

std::vector<double> Atom::electronDensityContinuous(const std::vector<double>& x) {
    std::vector<double> result(x.size(), 0.0);
    electronDensityContinuous(x.data(), result.data(), x.size());
    return result;
}

double Atom::electronDensityDiscrete(double x) {
    double result;
    result = electronDensity(x) - electronDensityContinuous(x);
    return result;
}

void Atom::electronDensityDiscrete(const double* x, double* y, std::size_t n) {
    // to do
}

std::vector<double> Atom::electronDensityDiscrete(const std::vector<double>& x) {
    std::vector<double> result(x.size(), 0.0);
    electronDensityDiscrete(x.data(), result.data(), x.size());
    return result;
}

}
}