#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <average-atom-toolkit/atom/thomas-fermi.h>
#include <average-atom-toolkit/atom/semiclassic.h>

namespace py = pybind11;

PYBIND11_MODULE(_pyaatk_atom, m) {

    py::class_<aatk::atom::Atom>(m, "Atom")
        .def(py::init([](
            double V, 
            double T, 
            double Z, 
            double tolerance,
            int    meshSize
        ) {
            auto atom = new aatk::atom::Atom(
                V, T, Z, tolerance, meshSize
            );
            return atom;
        }),
            py::arg("V")             = 1.0, 
            py::arg("T")             = 1.0,
            py::arg("Z")             = 1.0,
            py::arg("tolerance")     = 1.e-6,
            py::arg("meshSize")      = 1600
        )
        .def_property_readonly("volume", [](aatk::atom::Atom& atom) -> double {
            return atom.volume();
        })
        .def_property_readonly("temperature", [](aatk::atom::Atom& atom) -> double {
            return atom.temperature();
        })
        .def_property_readonly("Znucleus", [](aatk::atom::Atom& atom) -> double {
            return atom.Znucleus();
        })
        .def_property_readonly("chemicalPotential", [](aatk::atom::Atom& atom) -> double {
            return atom.chemicalPotential();
        })
        .def_property_readonly("ZfreeElectrons", [](aatk::atom::Atom& atom) -> double {
            return atom.ZfreeElectrons();
        })
        .def_property_readonly("radius", [](aatk::atom::Atom& atom) -> double {
            return atom.radius();
        })
        .def("U", [](aatk::atom::Atom& atom, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.U(x);
        })
        .def("U", [](aatk::atom::Atom& atom, py::array_t<double> x) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            bool correct_input = true;
            decltype(size) i = 0;
            while (i < size && correct_input) {
                correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
            }
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            auto  y = py::array_t<double>(size);
            auto cy = y.mutable_data();
            atom.U(cx, cy, size);
            return y;
        })
        .def("xU", [](aatk::atom::Atom& atom, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.xU(x);
        })
        .def("xU", [](aatk::atom::Atom& atom, py::array_t<double> x) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            bool correct_input = true;
            decltype(size) i = 0;
            while (i < size && correct_input) {
                correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
            }
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            auto  y = py::array_t<double>(size);
            auto cy = y.mutable_data();
            atom.xU(cx, cy, size);
            return y;
        })
        .def("x2dU", [](aatk::atom::Atom& atom, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.x2dU(x);
        })
        .def("x2dU", [](aatk::atom::Atom& atom, py::array_t<double> x) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            bool correct_input = true;
            decltype(size) i = 0;
            while (i < size && correct_input) {
                correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
            }
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            auto  y = py::array_t<double>(size);
            auto cy = y.mutable_data();
            atom.x2dU(cx, cy, size);
            return y;
        })
        .def("electronDensity", [](aatk::atom::Atom& atom, double x, double eb) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.electronDensity(x, eb);
        }, 
            py::arg("x"), 
            py::arg("eb") = 1.e+20
        )
        .def("electronDensity", [](aatk::atom::Atom& atom, py::array_t<double> x, double eb) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            bool correct_input = true;
            decltype(size) i = 0;
            while (i < size && correct_input) {
                correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
            }
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            auto  y = py::array_t<double>(size);
            auto cy = y.mutable_data();
            // write density to cy
            atom.electronDensity(cx, cy, size, eb);
            return y;
        },
            py::arg("x"), 
            py::arg("eb") = 1.e+20
        )
        .def("waveFunction", [](aatk::atom::Atom& atom, py::array_t<double> x, double e, double l) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            bool correct_input = true;
            decltype(size) i = 0;
            while (i < size && correct_input) {
                correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
            }
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            auto  y = py::array_t<double>(size);
            auto cy = y.mutable_data();
            // write density to cy
            atom.waveFunction(cx, cy, size, e, l);
            return y;
        })

        .def("outerRP", [](aatk::atom::Atom& atom, double e, double lambda) -> double {
            return atom.outerRP(e, lambda);
        })
        .def("innerRP", [](aatk::atom::Atom& atom, double e, double lambda) -> double {
            return atom.innerRP(e, lambda);
        })
        .def("action", [](aatk::atom::Atom& atom, double e, double lambda) -> double {
            return atom.action(e, lambda);
        })
        .def("energyLevel", [](aatk::atom::Atom& atom, int n, int l) -> double {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            if (l < 0 && l >= n) 
                throw std::runtime_error("Incorrect input: quantum number l should be between 0 and n - 1");
            return atom.energyLevel(n, l);
        })
        .def("electronStatesDiscrete", [](aatk::atom::Atom& atom, int n, int l) -> double {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            if (l < 0 && l >= n) 
                throw std::runtime_error("Incorrect input: quantum number l should be between 0 and n - 1");
            return atom.electronStatesDiscrete(n, l);
        })
        .def("electronStatesDiscrete", [](aatk::atom::Atom& atom, int n) -> double {
            return atom.electronStatesDiscrete(n);
        })
        .def("electronStatesDiscrete", [](aatk::atom::Atom& atom) -> double {
            return atom.electronStatesDiscrete();
        })
    ;

    py::class_<aatk::atom::ThomasFermiAtom, aatk::atom::Atom>(m, "ThomasFermiAtom")
        .def(py::init([](
            double V, 
            double T, 
            double Z, 
            double tolerance,
            int    meshSize
        ) {
            auto atom = new aatk::atom::ThomasFermiAtom(
                V, T, Z, tolerance, meshSize
            );
            return atom;
        }),
            py::arg("V")             = 1.0, 
            py::arg("T")             = 1.0,
            py::arg("Z")             = 1.0,
            py::arg("tolerance")     = 1.e-6,
            py::arg("meshSize")      = 1600
        )
    ;

    py::class_<aatk::atom::SemiclassicAtom, aatk::atom::Atom>(m, "SemiclassicAtom")
        .def(py::init([](
            double V, 
            double T, 
            double Z, 
            double tolerance,
            int    meshSize,
            int    nmax,
            bool   useContinuous
        ) {
            auto atom = new aatk::atom::SemiclassicAtom(
                V, T, Z, tolerance, meshSize, nmax, useContinuous
            );
            return atom;
        }),
            py::arg("V")             = 1.0, 
            py::arg("T")             = 1.0,
            py::arg("Z")             = 1.0,
            py::arg("tolerance")     = 1.e-6,
            py::arg("meshSize")      = 1600,
            py::arg("nmax")          = 20,
            py::arg("useContinuous") = true
        )
        .def("update", [](aatk::atom::SemiclassicAtom& atom, double mixing) -> void {
            atom.update(mixing);
        }, 
            py::arg("mixing") = 0.75
        )
        .def_property_readonly("discreteLevelsNumber", [](aatk::atom::SemiclassicAtom& atom) -> double {
            return atom.discreteLevelsNumber();
        })
        .def_property_readonly("boundaryEnergyValue", [](aatk::atom::SemiclassicAtom& atom) -> double {
            return atom.boundaryEnergyValue();
        })
        .def("electronDensityContinuous", [](aatk::atom::SemiclassicAtom& atom, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.electronDensityContinuous(x);
        })
        .def("electronDensityContinuous", [](aatk::atom::SemiclassicAtom& atom, py::array_t<double> x) -> py::array {
            if (x.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions: should be 1D array");
            }
            // get c-array representation
            const double* cx = x.data();
            std::size_t size = x.size();
            bool correct_input = true;
            decltype(size) i = 0;
            while (i < size && correct_input) {
                correct_input = (cx[i] >= 0.0 && cx[i] <= 1.0); ++i;
            }
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            auto  y = py::array_t<double>(size);
            auto cy = y.mutable_data();
            // write density to cy
            atom.electronDensityContinuous(cx, cy, size);
            return y;
        })
        .def("energyFull", [](aatk::atom::SemiclassicAtom& atom) -> double {
            return atom.energyFull();
        })

        .def("internalEnergy", [](aatk::atom::SemiclassicAtom& atom) -> double {
            return atom.internalEnergy();
        })

        .def("entropy", [](aatk::atom::SemiclassicAtom& atom) -> double {
            return atom.entropy();
        })
        .def("pressure", [](aatk::atom::SemiclassicAtom& atom) -> double {
            return atom.pressure();
        })
    ;
}