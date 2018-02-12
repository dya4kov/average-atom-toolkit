#pragma once
#include <array>

namespace numtools {
namespace ODE {

// ODE dimension to use with ODE solver
typedef unsigned Dimension;

// Array type to use with ODE solver
template<Dimension dim>
using Array = ::std::array<double, dim>;

}
}