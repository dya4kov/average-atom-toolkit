#include <boost/python/numpy.hpp>
#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace aatk {
namespace TF {

class ElectronStates {
public:
    double call_nl(int n, int l) {
        if (n < 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: quantum number n < 1");
            ::bpy::throw_error_already_set();
        }
        if (l < 0 && l >= n) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: quantum number l should be between 0 and n - 1");
            ::bpy::throw_error_already_set();
        }
        return N(n, l);
    }
    double call_n(int n) {
        if (n < 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: quantum number n < 1");
            ::bpy::throw_error_already_set();
        }
        return N(n);
    }
    double discrete() {
        return N.discrete();
    }
    double discrete_e(double e) {
        return N.discrete(e);
    }
    double continuous() {
        return N.continuous();
    }
    double continuous_e(double e) {
        return N.continuous(e);
    }
    double eBoundary() {
        return N.eBoundary();
    }
    void setV(double V) { N.setV(V); }
    void setT(double T) { N.setT(T); }
    void setZ(double Z) { N.setZ(Z); }
    void setVTZ(double V, double T, double Z) { N.setVTZ(V, T, Z); }
    void setNmax(int n) { N.setNmax(n); }
    void setTolerance(double eps) { N.setTolerance(eps); }

#ifdef ENABLE_MULTITHREADING
    void setThreadsLimit(int Nthreads) { N.setThreadsLimit(Nthreads); }
#endif

    ::bnp::ndarray discrete_en(::bnp::ndarray const & e) {
        if (e.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (e.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        auto ce   = reinterpret_cast<double*>(e.get_data());
        auto size = e.shape(0);
        auto cN   = N.discrete(ce, size);
        return ::bnp::from_data(
            cN,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(size),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    ::bnp::ndarray continuous_en(::bnp::ndarray const & e) {
        if (e.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (e.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        auto ce   = reinterpret_cast<double*>(e.get_data());
        auto size = e.shape(0);
        auto cN   = N.continuous(ce, size);
        return ::bnp::from_data(
            cN,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(size),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
private:
    ::aatk::TF::ElectronStates N;
};

}
}
}

BOOST_PYTHON_MODULE(_PyElectronStates) {
    bnp::initialize();

    bpy::class_<py::aatk::TF::ElectronStates>("ElectronStates")

        .def("__call__",     &py::aatk::TF::ElectronStates::call_nl)

        .def("__call__",     &py::aatk::TF::ElectronStates::call_n)

        .def("discrete",     &py::aatk::TF::ElectronStates::discrete)

        .def("discrete",     &py::aatk::TF::ElectronStates::discrete_e)

        .def("discrete",     &py::aatk::TF::ElectronStates::discrete_en)

        .def("continuous",   &py::aatk::TF::ElectronStates::continuous)

        .def("continuous",   &py::aatk::TF::ElectronStates::continuous_e)

        .def("continuous",   &py::aatk::TF::ElectronStates::continuous_en)

        .def("eBoundary",    &py::aatk::TF::ElectronStates::eBoundary)
        
        .def("setV",         &py::aatk::TF::ElectronStates::setV)

        .def("setT",         &py::aatk::TF::ElectronStates::setT)

        .def("setZ",         &py::aatk::TF::ElectronStates::setZ)

        .def("setVTZ",       &py::aatk::TF::ElectronStates::setVTZ)

        .def("setNmax",      &py::aatk::TF::ElectronStates::setNmax)

        .def("setTolerance", &py::aatk::TF::ElectronStates::setTolerance)
#ifdef ENABLE_MULTITHREADING
        .def("setThreadsLimit", &py::aatk::TF::ElectronStates::setThreadsLimit)
#endif
    ;

}