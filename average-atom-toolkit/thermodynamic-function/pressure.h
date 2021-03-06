#pragma once
#include <cstddef>
#include <vector>

namespace aatk {
namespace tfunc {

class Pressure : ThermodynamicFunction {
public:
	virtual double operator()(const double& V, const double& T);
	virtual std::vector<double> operator() (const std::vector<double>& V, const std::vector<double>& T);
	virtual double* operator()(const double* V, const double* T, const std::size_t& vsize, const std::size_t& tsize);

};

}
}