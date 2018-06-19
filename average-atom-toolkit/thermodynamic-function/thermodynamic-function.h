#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace aatk {
namespace tfunc {

class ThermodynamicFunction {
public:

    virtual double operator()(const double& V, const double& T);
    virtual std::vector<double> operator() (const std::vector<double>& V, const std::vector<double>& T);
    virtual double* operator()(const double* V, const double* T, const std::size_t& vsize, const std::size_t& tsize);

    virtual std::string name();

    void setZ(const double& Z);
    void setA(const double& A);
    void setThreadsLimit(const std::size_t& N);
    void setTolerance(const double& eps);

protected:

	double Z;
	double A;
    double tolerance;
    std::size_t threadsLimit;

};

}
}