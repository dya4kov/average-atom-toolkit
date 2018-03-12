#include <boost/python/numpy.hpp>
#include <average-atom-toolkit/thomas-fermi/atom/potential.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace aatk {
namespace TF {

class Potential {
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
        double* cx   = reinterpret_cast<double*>(x.get_data());
        auto    size = x.shape(0);
        // preprocess array to check if it is between 0 and 1
        bool correct_input = true;
        decltype(size) i = 0;
        while (i < size && correct_input) {
            correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
        }
        if (!correct_input) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: x should be between 0 and 1");
            ::bpy::throw_error_already_set();
        }
        auto cy = phi(cx, size);
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(size),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    ::bnp::ndarray dx_ndarray(::bnp::ndarray const & x) {
        if (x.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (x.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* cx   = reinterpret_cast<double*>(x.get_data());
        auto    size = x.shape(0);
        // preprocess array to check if it is between 0 and 1
        bool correct_input = true;
        decltype(size) i = 0;
        while (i < size && correct_input) {
            correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
        }
        if (!correct_input) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: x should be between 0 and 1");
            ::bpy::throw_error_already_set();
        }
        auto cy = phi.dx(cx, size);
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(size),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    double call_double(double x) {
        bool correct_input = (x >= 0.0 && x <= 1.0);
        if (!correct_input) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: x should be between 0 and 1");
            ::bpy::throw_error_already_set();
        }
        return phi(x);
    }
    double dx_double(double x) {
        bool correct_input = (x >= 0.0 && x <= 1.0);
        if (!correct_input) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: x should be between 0 and 1");
            ::bpy::throw_error_already_set();
        }
        return phi.dx(x);
    }
    void setV(double V) { phi.setV(V); }
    void setT(double T) { phi.setT(T); }
    void setZ(double Z) { phi.setZ(Z); }
    void setVTZ(double V, double T, double Z) { phi.setVTZ(V, T, Z); }
    void setTolerance(double eps) { phi.setTolerance(eps); }
private:
    ::aatk::TF::Potential phi;
};

}
}
}

BOOST_PYTHON_MODULE(_PyPotential) {
    bnp::initialize();

    bpy::class_<py::aatk::TF::Potential>("Potential")
        
        .def("__call__",     &py::aatk::TF::Potential::call_ndarray)

        .def("__call__",     &py::aatk::TF::Potential::call_double)

        .def("dx",           &py::aatk::TF::Potential::dx_ndarray)

        .def("dx",           &py::aatk::TF::Potential::dx_double)
        
        .def("setV",         &py::aatk::TF::Potential::setV)

        .def("setT",         &py::aatk::TF::Potential::setT)

        .def("setZ",         &py::aatk::TF::Potential::setZ)

        .def("setVTZ",       &py::aatk::TF::Potential::setVTZ)

        .def("setTolerance", &py::aatk::TF::Potential::setTolerance)

        ;

}