#include <cmath>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

using namespace numtools::specfunc::FD;

DMHalf::DMHalf() {
	a[0] = 1.71446374704454E+7;    b[0] = 9.67282587452899E+6;
	a[1] = 3.88148302324068E+7;    b[1] = 2.87386436731785E+7;
	a[2] = 3.16743385304962E+7;    b[2] = 3.26070130734158E+7;
	a[3] = 1.14587609192151E+7;    b[3] = 1.77657027846367E+7;
	a[4] = 1.83696370756153E+6;    b[4] = 4.81648022267831E+6;
	a[5] = 1.14980998186874E+5;    b[5] = 6.13709569333207E+5;
	a[6] = 1.98276889924768E+3;    b[6] = 3.13595854332114E+4;
	a[7] = 1.00000000000000E+0;    b[7] = 4.35061725080755E+2;

	c[ 0] = -4.46620341924942E-15; d[ 0] = -2.23310170962369E-15;
	c[ 1] = -1.58654991146236E-12; d[ 1] = -7.94193282071464E-13;
	c[ 2] = -4.44467627042232E-10; d[ 2] = -2.22564376956228E-10;
	c[ 3] = -6.84738791621745E-8;  d[ 3] = -3.43299431079845E-8;
	c[ 4] = -6.64932238528105E-6;  d[ 4] = -3.33919612678907E-6;
	c[ 5] = -3.69976170193942E-4;  d[ 5] = -1.86432212187088E-4;
	c[ 6] = -1.12295393687006E-2;  d[ 6] = -5.69764436880529E-3;
	c[ 7] = -1.60926102124442E-1;  d[ 7] = -8.34904593067194E-2;
	c[ 8] = -8.52408612877447E-1;  d[ 8] = -4.78770844009440E-1;
	c[ 9] = -7.45519953763928E-1;  d[ 9] = -4.99759250374148E-1;
	c[10] =  2.98435207466372E+0;  d[10] =  1.86795964993052E+0;
	c[11] =  1.00000000000000E+0;  d[11] =  4.16485970495288E-1;
}

double DMHalf::value(const double& x) {
	double  up = 0.0,  down = 0.0;
	double Dup = 0.0, Ddown = 0.0;
	double xpow[12] = {0.0};
	double t = 0.0;

	const int m1 = 8, m2 = 12;
	const int k1 = 8, k2 = 12;

	int i;
	xpow[0] = 1.0;
	if (x < 2.0) {
		t = std::exp(x);
		for (i = 1; i < m1; ++i) {
			xpow[i] = xpow[i - 1]*t;
		}
		for (i = 0; i < m1; ++i) {
			up += xpow[i]*a[i];
			down += xpow[i]*b[i];
		}
		for (i = 1; i < m1; ++i) {
			Dup += i*xpow[i - 1]*a[i];
			Ddown += i*xpow[i - 1]*b[i];
		}
		return t*(up/down + Dup*t/down - up*Ddown*t/(down*down));
	}
	else {
		t = std::sqrt(x);
		for (i = 1; i < m2; ++i) {
			xpow[i] = xpow[i - 1]/(x*x);
		}
		for (i = 0; i < m2; ++i) {
			up += xpow[i]*c[i];
			down += xpow[i]*d[i];
		}
		for (i = 1; i < m2; ++i) {
			Dup += i*xpow[i - 1]*c[i];
			Ddown += i*xpow[i - 1]*d[i];
		}
		return 0.5*up/down/t - 2.0/(x*x*t)*(Dup/down - up*Ddown/(down*down));
	}
}