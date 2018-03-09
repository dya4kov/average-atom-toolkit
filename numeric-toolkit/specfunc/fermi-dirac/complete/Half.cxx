#include <cmath>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

using namespace numtk::specfunc::FD;

Half::Half() {
	a[0] = 5.75834152995465E+6;    b[0] = 6.49759261942269E+6;
	a[1] = 1.30964880355883E+7;    b[1] = 1.70750501625775E+7;
	a[2] = 1.07608632249013E+7;    b[2] = 1.69288134856160E+7;
	a[3] = 3.93536421893014E+6;    b[3] = 7.95192647756086E+6;
	a[4] = 6.42493233715640E+5;    b[4] = 1.83167424554505E+6;
	a[5] = 4.16031909245777E+4;    b[5] = 1.95155948326832E+5;
	a[6] = 7.77238678539648E+2;    b[6] = 8.17922106644547E+3;
	a[7] = 1.00000000000000E+0;    b[7] = 9.02129136642157E+1;

	c[ 0] = 4.85378381173415E-14;  d[ 0] =	7.28067571760518E-14;
	c[ 1] = 1.64429113030738E-11;  d[ 1] =  2.45745452167585E-11;
	c[ 2] = 3.76794942277806E-9;   d[ 2] =  5.62152894375277E-9;
	c[ 3] = 4.69233883900644E-7;   d[ 3] =  6.96888634549649E-7;
	c[ 4] = 3.40679845803144E-5;   d[ 4] =  5.02360015186394E-5;
	c[ 5] = 1.32212995937796E-3;   d[ 5] =  1.92040136756592E-3;
	c[ 6] = 2.60768398973913E-2;   d[ 6] =  3.66887808002874E-2;
	c[ 7] = 2.48653216266227E-1;   d[ 7] =  3.24095226486468E-1;
	c[ 8] = 1.08037861921488E+0;   d[ 8] =  1.16434871200131E+0;
	c[ 9] = 1.91247528779676E+0;   d[ 9] =  1.34981244060549E+0;
	c[10] = 1.00000000000000E+0;   d[10] =  2.01311836975930E-1;
                                   d[11] = -2.14562434782759E-2;
}

double Half::value(const double& x) {
	double up = 0.0, down = 0.0;
	double xpow[12] = {0.0};
	double t = 0.0;

	const int m1 = 8, m2 = 11;
	const int k1 = 8, k2 = 12;

	int i;
	xpow[0] = 1.0;
	if (x < 2.0) {
		t = std::exp(x);
		for (i = 1; i < k1; ++i) {
			xpow[i] = xpow[i - 1]*t;
		}
		for (i = 0; i < m1; ++i) {
			up += xpow[i]*a[i];
		}
		for (i = 0; i < k1; ++i) {
			down += xpow[i]*b[i];
		}	
	}
	else {
		t = std::sqrt(x);
		for (i = 1; i < k2; ++i) {
			xpow[i] = xpow[i - 1]/(x*x);
		}
		for (i = 0; i < m2; ++i) {
			up += xpow[i]*c[i];			
		}
		for (i = 0; i < k2; ++i) {
			down += xpow[i]*d[i];
		}
		t = t*x;
	}	
	return t*up/down;
}