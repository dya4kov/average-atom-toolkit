#pragma once
#include <array>

namespace numtk {
namespace ODE {

// ODE dimension to use with ODE solver
typedef unsigned Dimension;

// Array type to use with ODE solver
template<Dimension dim>
using Array = ::std::array<double, dim>;

}
}