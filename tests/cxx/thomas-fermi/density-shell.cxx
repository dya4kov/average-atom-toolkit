#include <iostream>
#include <cmath>
#include <average-atom-toolkit/thomas-fermi/atom/shell/electron-density.h>

int main() {
	aatk::TF::shell::ElectronDensity rho;

	double V = 100.0;
	double T = 1.0;
	double Z = 78.0;

	rho.setVTZ(V, T, Z);

	int nx = 100;
	double* x = new double[nx];
	double xmin = 0.01; double xmax = 1.0;

    for (int ix = 0; ix < nx; ++ix)  {
        x[ix] = xmin + ix*(xmax - xmin)/(nx - 1);
        x[ix] *= x[ix];
    }

    auto result = rho(x, nx);

    double r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);

	for (int ix = 0; ix < nx; ++ix) {
		std::cout << x[ix] << ", " << 4.0*M_PI*r0*r0*x[ix]*x[ix]*result[ix] << std::endl;
	}

	delete[] result;

	return 0;
}