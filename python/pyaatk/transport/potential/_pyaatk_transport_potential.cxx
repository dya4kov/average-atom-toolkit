#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/transport/potential/debye.h>
#include <average-atom-toolkit/transport/potential/hard-sphere.h>
#include <average-atom-toolkit/transport/potential/polarization.h>
#include <average-atom-toolkit/transport/potential/square-well.h>
#include <average-atom-toolkit/transport/potential/tabular.h>

namespace py = pybind11;

using ::aatk::transport::potential::Base;
using ::aatk::transport::potential::Debye;
using ::aatk::transport::potential::HardSphere;
using ::aatk::transport::potential::Polarization;
using ::aatk::transport::potential::SquareWell;
using ::aatk::transport::potential::Tabular;

PYBIND11_MODULE(_pyaatk_transport_potential, m) {

	py::class_<Base>(m, "Base")
		.def("__call__", [](Base& potential, py::array_t<double> r) -> py::array {
    	    if (r.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            std::size_t size = r.size();
            auto cr = r.mutable_data();
            auto  V = py::array_t<double>(size);
            auto cV = V.mutable_data();
            potential(cr, cV, size);
            return V;
    	})
        .def("delta_eps", [](Base& potential, int l, double k) -> double  {
            return potential.delta_eps(l, k);
        });

	py::class_<Debye, Base>(m, "Debye")
		.def(py::init([](double Z, double eps) {
    	    auto potential = new Debye(Z, eps);
    	    return potential;
    	}), 
    	    py::arg("Z") = 1.0,
    	    py::arg("eps") = 1.e-6
    	);

    py::class_<HardSphere, Base>(m, "HardSphere")
		.def(py::init([](double eps) {
    	    auto potential = new HardSphere(eps);
    	    return potential;
    	}), 
    	    py::arg("eps") = 1.e-6
    	);

    py::class_<Polarization, Base>(m, "Polarization")
		.def(py::init([](double r0, double alpha, double eps) {
    	    auto potential = new Polarization(r0, alpha, eps);
    	    return potential;
    	}), 
    	    py::arg("r0")    = 1.0,
    	    py::arg("alpha") = 1.0,
    	    py::arg("eps") = 1.e-6
    	);

	py::class_<SquareWell, Base>(m, "SquareWell")
		.def(py::init([](double width, double depth, double eps) {
    	    auto potential = new SquareWell(width, depth, eps);
    	    return potential;
    	}), 
    	    py::arg("width") = 1.0,
    	    py::arg("depth") = 1.0,
    	    py::arg("eps") = 1.e-6
    	);

    py::class_<Tabular, Base>(m, "Tabular")
		.def(py::init([](py::array_t<double> r, py::array_t<double> rV, double eps) {
			if (r.size() != rV.size()) {
                throw std::runtime_error("r and rV arrays have different sizes");
            }
            std::size_t size = r.size();
            auto cr  = r.mutable_data();
            auto crV = rV.mutable_data();
    	    auto potential = new Tabular(cr, crV, size, eps);
    	    return potential;
    	}), 
    	    py::arg("r"),
    	    py::arg("rV"),
    	    py::arg("eps") = 1.e-6
    	);
}