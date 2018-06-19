#pragma once
#include <cstddef>
#include <vector>

#include <average-atom-toolkit/thermodynamic-function/free-energy.h>
#include <average-atom-toolkit/thermodynamic-function/chemical-potential.h>

namespace aatk {
namespace model {

class Model {
public:
	Model();
	virtual ::aatk::tfunc::FreeEnergy        FreeEnergy();
	virtual ::aatk::tfunc::ChemicalPotential ChemicalPotential();
};

}
}