#include <boost/python/numpy.hpp>
#include <average-atom-toolkit/thomas-fermi/thermodynamics/shell-chemical-potential.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace aatk {
namespace TF {
namespace shell {

class ChemicalPotential {
public:
    double call_double_double(double _V, double _T) { return M(_V,_T); };

    ::bnp::ndarray call_double_ndarray(double _V, ::bnp::ndarray _T) { 
        if (_T.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (_T.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* V   = new double[1]; V[0] = _V;
        double* T   = reinterpret_cast<double*>(_T.get_data());
        auto  tsize = _T.shape(0);
        decltype(tsize) vsize = 1;
        auto  data = M(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(tsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray call_ndarray_double(::bnp::ndarray _V, double _T) {
        if (_V.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (_V.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* V   = reinterpret_cast<double*>(_V.get_data());
        double* T   = new double[1]; T[0] = _T;
        auto  vsize = _V.shape(0);
        decltype(vsize) tsize = 1;
        auto  data = M(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray call_ndarray_ndarray(::bnp::ndarray _V, ::bnp::ndarray _T) { 
        if (_V.get_dtype() != ::bnp::dtype::get_builtin<double>() ||
            _T.get_dtype() != ::bnp::dtype::get_builtin<double>()  ) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (_V.get_nd() != 1 ||
            _T.get_nd() != 1  ) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* V  = reinterpret_cast<double*>(_V.get_data());
        double* T  = reinterpret_cast<double*>(_T.get_data());
        auto vsize = _V.shape(0);
        auto tsize = _T.shape(0);
        auto data  = M(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize, tsize),
            ::bpy::make_tuple(sizeof(double)*tsize, sizeof(double)),
            ::bpy::object()
        );
    };

    void setZ(double Z) { M.setZ(Z); }
    void setTolerance(double eps) { M.setTolerance(eps); }
    void setThreadsLimit(double Nthreads) { M.setThreadsLimit(Nthreads); }

private:
    ::aatk::TF::shell::ChemicalPotential M;
};

}
}
}
}

BOOST_PYTHON_MODULE(_PyChemicalPotentialShell) {
    bnp::initialize();

    bpy::class_<py::aatk::TF::shell::ChemicalPotential>("ChemicalPotential")
        
        .def("__call__",        &py::aatk::TF::shell::ChemicalPotential::call_double_double)
        .def("__call__",        &py::aatk::TF::shell::ChemicalPotential::call_double_ndarray)
        .def("__call__",        &py::aatk::TF::shell::ChemicalPotential::call_ndarray_double)
        .def("__call__",        &py::aatk::TF::shell::ChemicalPotential::call_ndarray_ndarray)

        .def("setZ",            &py::aatk::TF::shell::ChemicalPotential::setZ)
        .def("setTolerance",    &py::aatk::TF::shell::ChemicalPotential::setTolerance)
        .def("setThreadsLimit", &py::aatk::TF::shell::ChemicalPotential::setThreadsLimit)
    ;
}