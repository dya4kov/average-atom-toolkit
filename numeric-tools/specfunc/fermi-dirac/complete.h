#pragma once
/**
* Fermi-Dirac function approximations are taken from:
* H. M. Antia, "Rational function approximations for Fermi-Dirac Integrals",
* The Astrophysical Journal Supplement Series, 84:101-108, 1993
* The relative accuracy is up to the 1e-12
*/
namespace numtools {
namespace specfunc {

template<typename Func>
struct FermiDirac {
	FermiDirac() : f() {}
	double operator() (const double& x) {
		return f.value(x);
	}
private:
	Func f;
};

namespace FD {

struct DMHalf { 
	DMHalf();
	double value(const double& x);
private:
	double a[ 8]; 
	double b[ 8]; 
	double c[12]; 
	double d[12]; 
};
struct MHalf { 
	MHalf();
	double value(const double& x); 
private:
	double a[ 8];
	double b[ 8];
	double c[12];
	double d[12];
};
struct Half { 
	Half();
	double value(const double& x);
private:	
	double a[ 8];
	double b[ 8];
	double c[11];
	double d[12];
};
struct ThreeHalf { 
	ThreeHalf();
	double value(const double& x); 
private:
	double a[ 7];
	double b[ 8];
	double c[10];
	double d[11];
};

}
}
}