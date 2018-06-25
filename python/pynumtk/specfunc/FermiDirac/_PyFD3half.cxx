#include <boost/python/numpy.hpp>
#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace numtk {
namespace specfunc {
namespace FD {

class ThreeHalf {
public:
    ::bnp::ndarray call_ndarray(::bnp::ndarray const & x) {
        if (x.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (x.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        auto cx   = reinterpret_cast<double*>(x.get_data());
        auto size = x.shape(0);
        auto cy   = new double[size];
        for (decltype(size) i = 0; i < size; ++i) {
            cy[i] = FD3half.value(cx[i]);
        }
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(size),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    double call_double(double x) {
        return FD3half.value(x);
    }
private:
    ::numtk::specfunc::FD::ThreeHalf FD3half;
};

}
}
}
}

BOOST_PYTHON_MODULE(_PyFD3half) {
    bnp::initialize();

    bpy::class_<py::numtk::specfunc::FD::ThreeHalf>("ThreeHalf")
        
        .def("__call__", &py::numtk::specfunc::FD::ThreeHalf::call_ndarray)

        .def("__call__", &py::numtk::specfunc::FD::ThreeHalf::call_double)

        ;

}