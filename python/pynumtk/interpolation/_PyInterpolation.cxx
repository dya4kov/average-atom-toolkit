#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <numeric-toolkit/interpolation/linear.h>
#include <numeric-toolkit/interpolation/polynomial.h>
#include <numeric-toolkit/interpolation/spline.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyInterpolation, m) {

    py::class_<numtk::interpolation::Linear>(m, "Linear")
        .def(py::init([](py::array_t<double>& x, py::array_t<double>& y) {
        	if (x.ndim() != 1 || y.ndim() != 1) {
        		throw std::runtime_error("x and y should be 1D arrays");
        	}
        	if (x.shape(0) != y.shape(0)) {
        		throw std::runtime_error("x and y should have the same shape");
        	}
            auto interp = new numtk::interpolation::Linear(x.data(), y.data(), x.shape(0));
            return interp;
        }))

        .def("__call__", [](numtk::interpolation::Linear& interp, double x) -> double {
            return interp(x);
        })

        .def("__call__", [](numtk::interpolation::Linear& interp, py::array_t<double> x) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: x should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            auto y = py::array_t<double>(size);
            double* cy = y.mutable_data();
            for (decltype(size) i = 0; i < size; ++i) {
            	cy[i] = interp(cx[i]);
        	}
        	return y;
        })
    ;

    py::class_<numtk::interpolation::Polynomial>(m, "Polynomial")
        .def(py::init([](py::array_t<double>& x, py::array_t<double>& y, std::size_t order) {
        	if (x.ndim() != 1 || y.ndim() != 1) {
        		throw std::runtime_error("x and y should be 1D arrays");
        	}
        	if (x.shape(0) != y.shape(0)) {
        		throw std::runtime_error("x and y should have the same shape");
        	}
            auto interp = new numtk::interpolation::Polynomial(x.data(), y.data(), x.shape(0), order);
            return interp;
        }))

        .def("__call__", [](numtk::interpolation::Polynomial& interp, double x) -> double {
            return interp(x);
        })

        .def("__call__", [](numtk::interpolation::Polynomial& interp, py::array_t<double> x) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: x should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            auto y = py::array_t<double>(size);
            double* cy = y.mutable_data();
            for (decltype(size) i = 0; i < size; ++i) {
            	cy[i] = interp(cx[i]);
        	}
        	return y;
        })
    ;

    py::class_<numtk::interpolation::Spline>(m, "Spline")
        .def(py::init([](py::array_t<double>& x, py::array_t<double>& y) {
        	if (x.ndim() != 1 || y.ndim() != 1) {
        		throw std::runtime_error("x and y should be 1D arrays");
        	}
        	if (x.shape(0) != y.shape(0)) {
        		throw std::runtime_error("x and y should have the same shape");
        	}
            auto interp = new numtk::interpolation::Spline(x.data(), y.data(), x.shape(0));
            return interp;
        }))

        .def("__call__", [](numtk::interpolation::Spline& interp, double x) -> double {
            return interp(x);
        })

        .def("__call__", [](numtk::interpolation::Spline& interp, py::array_t<double> x) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: x should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            auto y = py::array_t<double>(size);
            double* cy = y.mutable_data();
            for (decltype(size) i = 0; i < size; ++i) {
            	cy[i] = interp(cx[i]);
        	}
        	return y;
        })
    ;
}