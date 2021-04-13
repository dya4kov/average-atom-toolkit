#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>


#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>

#include <average-atom-toolkit/atom/semiclassic.h>

namespace aatk {
namespace atom {

void SemiclassicAtom::electronDensity(const double* x, double* dens, std::size_t n, double eb) {
	for (std::size_t i = 0; i < n; ++i) {
		dens[i] = gsl_spline_eval(densSpline, std::sqrt(x[i]), densAcc);
	}
	return;
}

double SemiclassicAtom::electronDensity(double x) {
	double dens;
	electronDensity(&x, &dens, 1);
	return dens;
}

double SemiclassicAtom::electronDensityContinuous(double x) {
	numtk::specfunc::FermiDirac<numtk::specfunc::FD::Half>        FD_Half;
	numtk::specfunc::FermiDiracInc<numtk::specfunc::FDI::HalfInc> FD_Half_Inc;
	const double E_hartre = 27.21;

	const double E0 = boundaryEnergy;//energyLevel(nmax,nmax - 1);
	double V_r ;
	U(&x,&V_r,1);
	const double mu = M;
	const double factor = pow(2 * T,3.0 / 2.0) / (2 * pow(M_PI, 2) );
	const double y0 = (V_r + E0)/T;
	double result = 0.0;

	if (y0 <= 0){
		result =  FD_Half((V_r + mu)/T);
	}
	else{
		result = FD_Half((V_r + mu)/T) - FD_Half_Inc((V_r + mu)/T,y0); // fermi_dirac_1/2 + fermi_dirac_1/2_incomplete
	}

	return result * factor;
}

void SemiclassicAtom::electronDensityContinuous(const double* x, double* y, std::size_t n) {
	for (int i = 0; i < n; i++){
		y[i] = electronDensityContinuous(x[i]);
	}
}

std::vector<double> SemiclassicAtom::electronDensityContinuous(const std::vector<double>& x) {
	std::vector<double> result(x.size(), 0.0);
	electronDensityContinuous(x.data(), result.data(), x.size());
	return result;
}


}
}

