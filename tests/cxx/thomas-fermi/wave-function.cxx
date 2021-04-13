#include <iostream>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/shell/wave-function.h>

int main() {
	aatk::TF::shell::WaveFunction R;
	aatk::TF::EnergyLevel e;

	double V = 100.0;
	double T = 1.0;
	double Z = 78.0;

	e.setVTZ(V, T, Z);
	R.setVTZ(V, T, Z);

	double r0 = std::pow(3.0*V*Z / 4.0 / M_PI, 1.0 / 3.0);
	int n = 4; int l  = 0;

	double lambda = l + 0.5;
	double exact = M_PI*(n - l - 0.5);
	std::cout << "act = " << exact << std::endl;

	double lArg = 0.5*lambda*lambda / r0 / r0 * std::pow(Z, -2.0/3.0);
	double eArg = e(n,l)*std::pow(Z, -4.0/3.0);

	int nx = 100;
	double* x = new double[nx];
	double xmin = 0.01; double xmax = 1.0;

    for (int ix = 0; ix < nx; ++ix)  {
        x[ix] = xmin + ix*(xmax - xmin)/(nx - 1);
        x[ix] *= x[ix];
    }

    auto result = R(eArg, lArg, x, nx);

	for (int ix = 0; ix < nx; ++ix) {
		std::cout << x[ix] << ", " << result[ix] << std::endl;
	}

	return 0;
}