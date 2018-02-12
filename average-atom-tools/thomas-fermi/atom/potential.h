#pragma once
#include <vector>
#include <average-atom-tools/thomas-fermi/atom/table.h>

namespace AATools {
namespace TF {

class Potential {
public:
    Potential();
    Potential(const Potential& potential);
    Potential& operator=(const Potential& potential);
    
    double mu();

    double operator()(const double& _x);
    double dx(const double& _x);
    
    std::vector<double>& operator()(const std::vector<double>& _xarray);
    std::vector<double>& dx(const std::vector<double>& _xarray);

    double* operator()(double* _xarray, size_t n);
    double* dx(double* _xarray, size_t n);

    void setV(const double& V);
    void setT(const double& T);
    void setZ(const double& Z);

    void setTolerance(const double& eps);

private:

	double V1, T1, mu1;
    double VZ, TZ, muZ;
	double tolerance;
	bool   needUpdate;

    const double bestTolerance;

    Table table;
};

}
}