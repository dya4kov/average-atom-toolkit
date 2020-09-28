#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <average-atom-toolkit/semiclassic/atom.h>
#ifdef ENABLE_MULTITHREADING
#include <average-atom-toolkit/multithreading/thread-pool.h>
#endif

namespace py = pybind11;

PYBIND11_MODULE(_aatk_semiclassic_atom, m) {

    py::class_<aatk::semiclassic::Atom>(m, "Atom")
        .def(py::init([](
            double V, 
            double T, 
            double Z, 
            int    nmax,
            double tol
     #ifdef ENABLE_MULTITHREADING
            ,::aatk::multithreading::ThreadPool& pool
     #endif
        ) {
            auto atom = new aatk::semiclassic::Atom(V, T, Z, nmax, tol
     #ifdef ENABLE_MULTITHREADING
                ,pool
     #endif 
            );
            return atom;
        }),
            py::arg("V")         = 1.0, 
            py::arg("T")         = 1.0,
            py::arg("Z")         = 1.0,
            py::arg("nmax")      = 20,
            py::arg("tolerance") = 1.e-6
     #ifdef ENABLE_MULTITHREADING
            ,py::arg("threads") = ::aatk::multithreading::dummy_pool
     #endif
        )

        .def("reset", [](aatk::semiclassic::Atom& atom,
            double       V, 
            double       T, 
            double       Z, 
            int          nmax
        ) -> void {
            atom.reset();
        }, 
            py::arg("V")    = -1.0, 
            py::arg("T")    = -1.0,
            py::arg("Z")    = -1.0,
            py::arg("nmax") = -1
        )
        .def("update", [](aatk::semiclassic::Atom& atom, double mixing) -> void {
            atom.update(mixing);
        }, 
            py::arg("mixing") = 0.25
        )
        .def("update", [](aatk::semiclassic::Atom& atom, py::array_t<double> x, double mixing) -> void {
            atom.update(x.data(), x.size(), mixing);
        },
            py::arg("mesh"),
            py::arg("mixing") = 0.25
        )
        .def_property_readonly("V", [](aatk::semiclassic::Atom& atom) -> double {
            return atom.V();
        })
        .def_property_readonly("T", [](aatk::semiclassic::Atom& atom) -> double {
            return atom.T();
        })
        .def_property_readonly("Z", [](aatk::semiclassic::Atom& atom) -> double {
            return atom.Z();
        })
        .def_property_readonly("M", [](aatk::semiclassic::Atom& atom) -> double {
            return atom.M();
        })

        // .def("setTolerance", [](aatk::semiclassic::Atom& atom, double tol) {
        //     atom.setTolerance(tol);
        // })

        .def("energyLevel", [](aatk::semiclassic::Atom& atom, int n, int l) -> double {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            if (l < 0 && l >= n) 
                throw std::runtime_error("Incorrect input: quantum number l should be between 0 and n - 1");
            return atom.energyLevel(n, l);
        })

        // .def("energyLevel", [](aatk::semiclassic::Atom& atom, int n) -> py::array {
        //     if (n < 1) 
        //         throw std::runtime_error("Incorrect input: quantum number n < 1");
        //     auto pyen = py::array_t<double>(n);
        //     auto pydata = pyen.mutable_data();
        //     auto en = atom.energyLevel(n);
        //     for (int l = 0; l < n; ++l) pydata[l] = en[l];
        //     return pyen;
        // })

        .def("potential", [](aatk::semiclassic::Atom& atom, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.potential(x);
        })

        .def("potential", [](aatk::semiclassic::Atom& atom, py::array_t<double> x) -> py::array {
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
            atom.potential(cx, cy, size);
            return y;
        })

        .def("waveFunction", [](aatk::semiclassic::Atom& atom, double e, double l, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.waveFunction(e, l, x);
        })
        .def("waveFunction", [](aatk::semiclassic::Atom& atom, double e, double l, py::array_t<double> x) -> py::array {
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

            atom.waveFunction(e, l, cx, cy, size);
            return y;
        })

        .def("electronDensity", [](aatk::semiclassic::Atom& atom, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return atom.electronDensity(x);
        })
        .def("electronDensity", [](aatk::semiclassic::Atom& atom, py::array_t<double> x) -> py::array {
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

            atom.electronDensity(cx, cy, size);
            return y;
        })

        .def("electronStates", [](aatk::semiclassic::Atom& atom, int n, int l) -> double {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            if (l < 0 && l >= n) 
                throw std::runtime_error("Incorrect input: quantum number l should be between 0 and n - 1");
            return atom.electronStates(n, l);
        })
        .def("electronStates", [](aatk::semiclassic::Atom& atom, int n) -> double {
            // if (n < 1) 
            //     throw std::runtime_error("Incorrect input: quantum number n < 1");
            // auto pyN = py::array_t<double>(n);
            // auto pydata = pyN.mutable_data();
            // auto N = atom.electronStates(n);
            // for (int l = 0; l < n; ++l) pydata[l] = N[l];
            // return pyN;
            return atom.electronStates(n);
        })
        .def("electronStates", [](aatk::semiclassic::Atom& atom) -> double {
            return atom.electronStates();
        })

        .def("innerRP", [](aatk::semiclassic::Atom& atom, double e, double l) -> double {
            return atom.innerRP(e, l)[0];
        })
        .def("outerRP", [](aatk::semiclassic::Atom& atom, double e, double l) -> double {
            return atom.outerRP(e, l)[0];
        })
    ;
}