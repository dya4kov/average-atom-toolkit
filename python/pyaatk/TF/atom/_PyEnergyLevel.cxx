#include <boost/python/numpy.hpp>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace aatk {
namespace TF {

class EnergyLevel {
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
        return e(n, l);
    }
    bnp::ndarray getitem_n(int n) {
        if (n < 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect input: quantum number n < 1");
            ::bpy::throw_error_already_set();
        }
        double* cen = new double[n];
        auto     en = e[n];
        for (int l = 0; l < n; ++l) cen[l] = en[l];
        return ::bnp::from_data(
            cen,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(n),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    void setThreadsLimit(int Nthreads) { e.setThreadsLimit(std::max(1, Nthreads)); }
    void prepareLevelsBelow(int nMax) { e.prepareLevelsBelow(nMax); }
    void setV(double V) { e.setV(V); }
    void setT(double T) { e.setT(T); }
    void setZ(double Z) { e.setZ(Z); }
    void setVTZ(double V, double T, double Z) { e.setVTZ(V, T, Z); }
    void setTolerance(double eps) { e.setTolerance(eps); }
private:
    ::aatk::TF::EnergyLevel e;
};

}
}
}

BOOST_PYTHON_MODULE(_PyEnergyLevel) {
    bnp::initialize();

    bpy::class_<py::aatk::TF::EnergyLevel>("EnergyLevel")

        .def("__call__",           &py::aatk::TF::EnergyLevel::call_nl)

        .def("__getitem__",        &py::aatk::TF::EnergyLevel::getitem_n)

        .def("prepareLevelsBelow", &py::aatk::TF::EnergyLevel::prepareLevelsBelow)

        .def("setV",               &py::aatk::TF::EnergyLevel::setV)

        .def("setT",               &py::aatk::TF::EnergyLevel::setT)

        .def("setZ",               &py::aatk::TF::EnergyLevel::setZ)

        .def("setVTZ",             &py::aatk::TF::EnergyLevel::setVTZ)

        .def("setTolerance",       &py::aatk::TF::EnergyLevel::setTolerance)

    ;
}