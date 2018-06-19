#include <boost/python/numpy.hpp>
#include <average-atom-toolkit/thomas-fermi/atom/electron-density.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace aatk {
namespace TF {
namespace qx {

class ElectronDensity {
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
        auto cy = rho(cx, size);
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
        return rho(x);
    }
    void setV(double V) { rho.setV(V); }
    void setT(double T) { rho.setT(T); }
    void setZ(double Z) { rho.setZ(Z); }
    void setVTZ(double V, double T, double Z) { rho.setVTZ(V, T, Z); }
    void setTolerance(double eps) { rho.setTolerance(eps); }
    void setBoundary(double eb)   { rho.setBoundary(eb); }
private:
    ::aatk::TF::ElectronDensity rho;
};

}
}
}

BOOST_PYTHON_MODULE(_PyElectronDensity) {
    bnp::initialize();

    bpy::class_<py::aatk::TF::ElectronDensity>("ElectronDensity")
        
        .def("__call__",        &py::aatk::TF::ElectronDensity::call_ndarray)

        .def("__call__",        &py::aatk::TF::ElectronDensity::call_double)

        .def("setV",            &py::aatk::TF::ElectronDensity::setV)

        .def("setT",            &py::aatk::TF::ElectronDensity::setT)

        .def("setZ",            &py::aatk::TF::ElectronDensity::setZ)

        .def("setVTZ",          &py::aatk::TF::ElectronDensity::setVTZ)

        .def("setTolerance",    &py::aatk::TF::ElectronDensity::setTolerance)

        ;

}