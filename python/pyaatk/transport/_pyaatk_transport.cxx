#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/transport/potential/base.h>
#include <average-atom-toolkit/transport/cross-section.h>
#include <average-atom-toolkit/transport/Kfunction.h>

namespace py = pybind11;

using Potenial = ::aatk::transport::potential::Base;
using ::aatk::transport::CrossSection;
using ::aatk::transport::Kfunction;

PYBIND11_MODULE(_pyaatk_transport, m) {

	py::class_<CrossSection>(m, "CrossSection")
        .def(py::init([](Potenial& V, int nterms, double rmax) {
            auto cs = new CrossSection(V, nterms, rmax);
            return cs;
        }), 
            py::arg("potential"),
            py::arg("nterms") = 30,
            py::arg("rmax") = 10000.
        )
		.def("__call__", [](CrossSection& cs, py::array_t<double> k) -> py::array {
    	    if (k.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            std::size_t size = k.size();
            auto ck = k.mutable_data();
            auto sigma = py::array_t<double>(size);
            auto csigma = sigma.mutable_data();
            cs(ck, csigma, size);
            return sigma;
    	})
        .def("__call__", [](CrossSection& cs, double k) -> double {
            return cs(k);
        });

    py::class_<Kfunction>(m, "Kfunction")
        .def(py::init([](double n, double T, double M) {
            auto K = new Kfunction(n, T, M);
            return K;
        }), 
            py::arg("n"),
            py::arg("T"),
            py::arg("M")
        )
        .def("__call__", [](Kfunction& K, py::array_t<double> e, py::array_t<double> tau) -> double {
            if (e.size() != tau.size()) {
                throw std::runtime_error("arrays should be of the same shape");
            }
            // get c-array representation
            std::size_t size = e.size();
            auto ce = e.mutable_data();
            auto ctau = tau.mutable_data();
            return K(ce, ctau, size);
        });
}