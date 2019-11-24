#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyFermiDirac, m) {
	py::class_<numtk::specfunc::FD::ThreeHalf>(m, "ThreeHalf")
        .def(py::init([]() {
            auto FD = new numtk::specfunc::FD::ThreeHalf();
            return FD;
        }))
        .def("__call__", [](numtk::specfunc::FD::ThreeHalf& FD, double x) -> double {
            return FD.value(x);
        })
        .def("__call__", [](numtk::specfunc::FD::ThreeHalf& FD, py::array_t<double> x) -> py::array {
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            auto y = py::array_t<double>(x.request().shape);
            double* cy = y.mutable_data();
            for (std::size_t i = 0; i < size; ++i) {
            	cy[i] = FD.value(cx[i]);
        	}
        	return y;
        })
    ;

    py::class_<numtk::specfunc::FD::Half>(m, "Half")
        .def(py::init([]() {
            auto FD = new numtk::specfunc::FD::Half();
            return FD;
        }))
        .def("__call__", [](numtk::specfunc::FD::Half& FD, double x) -> double {
            return FD.value(x);
        })
        .def("__call__", [](numtk::specfunc::FD::Half& FD, py::array_t<double> x) -> py::array {
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            auto y = py::array_t<double>(x.request().shape);
            double* cy = y.mutable_data();
            for (std::size_t i = 0; i < size; ++i) {
            	cy[i] = FD.value(cx[i]);
        	}
        	return y;
        })
    ;

    py::class_<numtk::specfunc::FD::MHalf>(m, "MHalf")
        .def(py::init([]() {
            auto FD = new numtk::specfunc::FD::MHalf();
            return FD;
        }))
        .def("__call__", [](numtk::specfunc::FD::MHalf& FD, double x) -> double {
            return FD.value(x);
        })
        .def("__call__", [](numtk::specfunc::FD::MHalf& FD, py::array_t<double> x) -> py::array {
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            auto y = py::array_t<double>(x.request().shape);
            double* cy = y.mutable_data();
            for (std::size_t i = 0; i < size; ++i) {
            	cy[i] = FD.value(cx[i]);
        	}
        	return y;
        })
    ;

    py::class_<numtk::specfunc::FD::DMHalf>(m, "DMHalf")
        .def(py::init([]() {
            auto FD = new numtk::specfunc::FD::DMHalf();
            return FD;
        }))
        .def("__call__", [](numtk::specfunc::FD::DMHalf& FD, double x) -> double {
            return FD.value(x);
        })
        .def("__call__", [](numtk::specfunc::FD::DMHalf& FD, py::array_t<double> x) -> py::array {
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            auto y = py::array_t<double>(x.request().shape);
            double* cy = y.mutable_data();
            for (std::size_t i = 0; i < size; ++i) {
            	cy[i] = FD.value(cx[i]);
        	}
        	return y;
        })
    ;
}