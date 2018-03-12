#include <boost/python/numpy.hpp>
#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>

namespace bpy = boost::python;
namespace bnp = boost::python::numpy;

namespace py {
namespace aatk {
namespace TF {

class RotatePoints {
public:
    ::bnp::ndarray inner_ndarray_ndarray(::bnp::ndarray const & e, ::bnp::ndarray const & l) {
        if (e.get_dtype() != ::bnp::dtype::get_builtin<double>() ||
            l.get_dtype() != ::bnp::dtype::get_builtin<double>() ) 
        {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (e.get_nd() != 1 || l.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* ce = reinterpret_cast<double*>(e.get_data());
        double* cl = reinterpret_cast<double*>(l.get_data());
        auto    esize = e.shape(0);
        auto    lsize = l.shape(0);
        double* cy = new double[esize*lsize];
        for (decltype(esize) ie = 0; ie < esize; ++ie) {
            for (decltype(lsize) il = 0; il < lsize; ++il) {
                cy[ie*lsize + il] = RP.inner(ce[ie], cl[il]);
            }
        }
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(esize, lsize),
            ::bpy::make_tuple(sizeof(double)*lsize, sizeof(double)),
            ::bpy::object()
        );
    }
    ::bnp::ndarray outer_ndarray_ndarray(::bnp::ndarray const & e, ::bnp::ndarray const & l) {
        if (e.get_dtype() != ::bnp::dtype::get_builtin<double>() ||
            l.get_dtype() != ::bnp::dtype::get_builtin<double>() ) 
        {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (e.get_nd() != 1 || l.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* ce = reinterpret_cast<double*>(e.get_data());
        double* cl = reinterpret_cast<double*>(l.get_data());
        auto    esize = e.shape(0);
        auto    lsize = l.shape(0);
        double* cy = new double[esize*lsize];
        for (decltype(esize) ie = 0; ie < esize; ++ie) {
            for (decltype(lsize) il = 0; il < lsize; ++il) {
                cy[ie*lsize + il] = RP.outer(ce[ie], cl[il]);
            }
        }
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(esize, lsize),
            ::bpy::make_tuple(sizeof(double)*lsize, sizeof(double)),
            ::bpy::object()
        );
    }
    ::bnp::ndarray inner_ndarray_double(::bnp::ndarray const & e, double l) {
        if (e.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (e.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* ce = reinterpret_cast<double*>(e.get_data());
        auto    esize = e.shape(0);
        double* cy = new double[esize];
        for (decltype(esize) ie = 0; ie < esize; ++ie) {
            cy[ie] = RP.inner(ce[ie], l);
        }
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(esize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    ::bnp::ndarray outer_ndarray_double(::bnp::ndarray const & e, double l) {
        if (e.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (e.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* ce = reinterpret_cast<double*>(e.get_data());
        auto    esize = e.shape(0);
        double* cy = new double[esize];
        for (decltype(esize) ie = 0; ie < esize; ++ie) {
            cy[ie] = RP.outer(ce[ie], l);
        }
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(esize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    ::bnp::ndarray inner_double_ndarray(double e, ::bnp::ndarray const & l) {
        if (l.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (l.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* cl = reinterpret_cast<double*>(l.get_data());
        auto    lsize = l.shape(0);
        double* cy = new double[lsize];
        for (decltype(lsize) il = 0; il < lsize; ++il) {
            cy[il] = RP.inner(e, cl[il]);
        }
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(lsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    ::bnp::ndarray outer_double_ndarray(double e, ::bnp::ndarray const & l) {
        if (l.get_dtype() != ::bnp::dtype::get_builtin<double>()) {
            PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
            ::bpy::throw_error_already_set();
        }
        if (l.get_nd() != 1) {
            PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
            ::bpy::throw_error_already_set();
        }
        // get c-array representation
        double* cl = reinterpret_cast<double*>(l.get_data());
        auto    lsize = l.shape(0);
        double* cy = new double[lsize];
        for (decltype(lsize) il = 0; il < lsize; ++il) {
            cy[il] = RP.outer(e, cl[il]);
        }
        return ::bnp::from_data(
            cy,
            ::bnp::dtype::get_builtin<double>(),
            ::bpy::make_tuple(lsize),
            ::bpy::make_tuple(sizeof(double)),
            ::bpy::object()
        );
    }
    double inner_double_double(double e, double l) {
        return RP.inner(e, l);
    }
    double outer_double_double(double e, double l) {
        return RP.outer(e, l);
    }
    void setV(double V) { RP.setV(V); }
    void setT(double T) { RP.setT(T); }
    void setZ(double Z) { RP.setZ(Z); }
    void setVTZ(double V, double T, double Z) { RP.setVTZ(V, T, Z); }
    void setTolerance(double eps) { RP.setTolerance(eps); }
private:
    ::aatk::TF::RotatePoints RP;
};

}
}
}

BOOST_PYTHON_MODULE(_PyRotatePoints) {
    bnp::initialize();

    bpy::class_<py::aatk::TF::RotatePoints>("RotatePoints")
        
        .def("inner",        &py::aatk::TF::RotatePoints::inner_ndarray_ndarray)

        .def("outer",        &py::aatk::TF::RotatePoints::outer_ndarray_ndarray)

        .def("inner",        &py::aatk::TF::RotatePoints::inner_double_ndarray)

        .def("outer",        &py::aatk::TF::RotatePoints::outer_double_ndarray)

        .def("inner",        &py::aatk::TF::RotatePoints::inner_ndarray_double)

        .def("outer",        &py::aatk::TF::RotatePoints::outer_ndarray_double)

        .def("inner",        &py::aatk::TF::RotatePoints::inner_double_double)

        .def("outer",        &py::aatk::TF::RotatePoints::outer_double_double)
        
        .def("setV",         &py::aatk::TF::RotatePoints::setV)

        .def("setT",         &py::aatk::TF::RotatePoints::setT)

        .def("setZ",         &py::aatk::TF::RotatePoints::setZ)

        .def("setVTZ",       &py::aatk::TF::RotatePoints::setVTZ)

        .def("setTolerance", &py::aatk::TF::RotatePoints::setTolerance)

    ;

}