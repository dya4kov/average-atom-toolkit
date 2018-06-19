#pragma once
#include <cstddef>
#include <vector>

#include <average-atom-toolkit/thermodynamic-function/free-energy.h>
#include <average-atom-toolkit/thermodynamic-function/chemical-potential.h>

#include <average-atom-toolkit/thomas-fermi/eos/free-energy.h>

namespace aatk {
namespace model {

class ThomasFermi : Model {
	::aatk::tfunc::FreeEnergy FreeEnergy() {
		return ::aatk::TF::FreeEnergy();
	}
	::aatk::tfunc::ChemicalPotential ChemicalPotential() {
		return ::aatk::TF::ChemicalPotential();
	}
};

}
}