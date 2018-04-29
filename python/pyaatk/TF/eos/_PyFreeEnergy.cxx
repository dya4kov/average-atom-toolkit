#include <boost/python/numpy.hpp>
#include <average-atom-toolkit/thomas-fermi/eos/free-energy.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace aatk {
namespace TF {

class FreeEnergy {
public:
    double call_double_double(double _V, double _T) { return F(_V,_T); };

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
        auto  data = F(V, T, vsize, tsize);
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
        auto  data = F(V, T, vsize, tsize);
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
        auto data  = F(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize, tsize),
            ::bpy::make_tuple(sizeof(double)*tsize, sizeof(double)),
            ::bpy::object()
        );
    };

    double DV_double_double(double _V, double _T) { return F.DV(_V,_T); };

    ::bnp::ndarray DV_double_ndarray(double _V, ::bnp::ndarray _T) { 
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
        auto  data = F.DV(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(tsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray DV_ndarray_double(::bnp::ndarray _V, double _T) {
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
        auto  data = F.DV(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray DV_ndarray_ndarray(::bnp::ndarray _V, ::bnp::ndarray _T) { 
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
        auto data  = F.DV(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize, tsize),
            ::bpy::make_tuple(sizeof(double)*tsize, sizeof(double)),
            ::bpy::object()
        );
    };

    double DT_double_double(double _V, double _T) { return F.DT(_V,_T); };

    ::bnp::ndarray DT_double_ndarray(double _V, ::bnp::ndarray _T) { 
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
        auto  data = F.DT(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(tsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray DT_ndarray_double(::bnp::ndarray _V, double _T) {
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
        auto  data = F.DT(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray DT_ndarray_ndarray(::bnp::ndarray _V, ::bnp::ndarray _T) { 
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
        auto data  = F.DT(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize, tsize),
            ::bpy::make_tuple(sizeof(double)*tsize, sizeof(double)),
            ::bpy::object()
        );
    };

    double D2V_double_double(double _V, double _T) { return F.D2V(_V,_T); };

    ::bnp::ndarray D2V_double_ndarray(double _V, ::bnp::ndarray _T) { 
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
        auto  data = F.D2V(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(tsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray D2V_ndarray_double(::bnp::ndarray _V, double _T) {
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
        auto  data = F.D2V(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray D2V_ndarray_ndarray(::bnp::ndarray _V, ::bnp::ndarray _T) { 
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
        auto data  = F.D2V(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize, tsize),
            ::bpy::make_tuple(sizeof(double)*tsize, sizeof(double)),
            ::bpy::object()
        );
    };

    double DVT_double_double(double _V, double _T) { return F.DVT(_V,_T); };

    ::bnp::ndarray DVT_double_ndarray(double _V, ::bnp::ndarray _T) { 
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
        auto  data = F.DVT(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(tsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray DVT_ndarray_double(::bnp::ndarray _V, double _T) {
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
        auto  data = F.DVT(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray DVT_ndarray_ndarray(::bnp::ndarray _V, ::bnp::ndarray _T) { 
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
        auto data  = F.DVT(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize, tsize),
            ::bpy::make_tuple(sizeof(double)*tsize, sizeof(double)),
            ::bpy::object()
        );
    };

    double D2T_double_double(double _V, double _T) { return F.D2T(_V,_T); };

    ::bnp::ndarray D2T_double_ndarray(double _V, ::bnp::ndarray _T) { 
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
        auto  data = F.D2T(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(tsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray D2T_ndarray_double(::bnp::ndarray _V, double _T) {
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
        auto  data = F.D2T(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    };
    ::bnp::ndarray D2T_ndarray_ndarray(::bnp::ndarray _V, ::bnp::ndarray _T) { 
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
        auto data  = F.D2T(V, T, vsize, tsize);
        return ::bnp::from_data(
            data,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(vsize, tsize),
            ::bpy::make_tuple(sizeof(double)*tsize, sizeof(double)),
            ::bpy::object()
        );
    };

    void setZ(double Z) { F.setZ(Z); }
    void setTolerance(double eps) { F.setTolerance(eps); }
    void setThreadsLimit(double Nthreads) { F.setThreadsLimit(Nthreads); }

private:
    ::aatk::TF::FreeEnergy F;
};
}
}
}

BOOST_PYTHON_MODULE(_PyFreeEnergy) {
    bnp::initialize();

    bpy::class_<py::aatk::TF::FreeEnergy>("FreeEnergy")
        
        .def("__call__",        &py::aatk::TF::FreeEnergy::call_double_double)
        .def("__call__",        &py::aatk::TF::FreeEnergy::call_double_ndarray)
        .def("__call__",        &py::aatk::TF::FreeEnergy::call_ndarray_double)
        .def("__call__",        &py::aatk::TF::FreeEnergy::call_ndarray_ndarray)

        .def("DV",              &py::aatk::TF::FreeEnergy::DV_double_double)
        .def("DV",              &py::aatk::TF::FreeEnergy::DV_double_ndarray)
        .def("DV",              &py::aatk::TF::FreeEnergy::DV_ndarray_double)
        .def("DV",              &py::aatk::TF::FreeEnergy::DV_ndarray_ndarray)

        .def("DT",              &py::aatk::TF::FreeEnergy::DT_double_double)
        .def("DT",              &py::aatk::TF::FreeEnergy::DT_double_ndarray)
        .def("DT",              &py::aatk::TF::FreeEnergy::DT_ndarray_double)
        .def("DT",              &py::aatk::TF::FreeEnergy::DT_ndarray_ndarray)

        .def("D2V",             &py::aatk::TF::FreeEnergy::D2V_double_double)
        .def("D2V",             &py::aatk::TF::FreeEnergy::D2V_double_ndarray)
        .def("D2V",             &py::aatk::TF::FreeEnergy::D2V_ndarray_double)
        .def("D2V",             &py::aatk::TF::FreeEnergy::D2V_ndarray_ndarray)

        .def("DVT",             &py::aatk::TF::FreeEnergy::DVT_double_double)
        .def("DVT",             &py::aatk::TF::FreeEnergy::DVT_double_ndarray)
        .def("DVT",             &py::aatk::TF::FreeEnergy::DVT_ndarray_double)
        .def("DVT",             &py::aatk::TF::FreeEnergy::DVT_ndarray_ndarray)

        .def("D2T",             &py::aatk::TF::FreeEnergy::D2T_double_double)
        .def("D2T",             &py::aatk::TF::FreeEnergy::D2T_double_ndarray)
        .def("D2T",             &py::aatk::TF::FreeEnergy::D2T_ndarray_double)
        .def("D2T",             &py::aatk::TF::FreeEnergy::D2T_ndarray_ndarray)

        .def("setZ",            &py::aatk::TF::FreeEnergy::setZ)
        .def("setTolerance",    &py::aatk::TF::FreeEnergy::setTolerance)
        .def("setThreadsLimit", &py::aatk::TF::FreeEnergy::setThreadsLimit)
    ;
}