#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/thomas-fermi/atom/electron-density.h>
#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>
#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>
#include <average-atom-toolkit/thomas-fermi/atom/potential.h>
#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>

namespace py = pybind11;

PYBIND11_MODULE(_pyaatk_thomas_fermi_atom, m) {

    auto& api = py::detail::npy_api::get();

    py::class_<aatk::TF::ElectronDensity>(m, "ElectronDensity")
        .def(py::init([](){
            auto rho = new aatk::TF::ElectronDensity();
            return rho;
        }))

        .def("__call__", [](aatk::TF::ElectronDensity& rho, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return rho(x);
        })

        .def("__call__", [api](aatk::TF::ElectronDensity& rho, py::array_t<double> x) -> py::array {
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
            auto cy = rho(cx, size);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          x.ndim(), 
                    (Py_intptr_t*) x.shape(), 
                    (Py_intptr_t*) x.strides(),
                    (void *)       cy, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("setV", [](aatk::TF::ElectronDensity& rho, double V) {
            rho.setV(V);
        })
        .def("setT", [](aatk::TF::ElectronDensity& rho, double T) {
            rho.setT(T);
        })
        .def("setZ", [](aatk::TF::ElectronDensity& rho, double Z) {
            rho.setZ(Z);
        })
        .def("setVTZ", [](aatk::TF::ElectronDensity& rho, double V, double T, double Z) {
            rho.setVTZ(V,T,Z);
        })
        .def("setTolerance", [](aatk::TF::ElectronDensity& rho, double tol){
            rho.setTolerance(tol);
        })
        .def("setBoundary", [](aatk::TF::ElectronDensity& rho, double eb){
            rho.setBoundary(eb);
        })
    ;

    py::class_<aatk::TF::ElectronStates>(m, "ElectronStates")
        .def(py::init([](){
            auto N = new aatk::TF::ElectronStates();
            return N;
        }))

        .def("__call__", [](aatk::TF::ElectronStates& N, int n, int l) -> double {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            if (l < 0 && l >= n) 
                throw std::runtime_error("Incorrect input: quantum number l should be between 0 and n - 1");
            return N(n, l);
        })
        .def("__call__", [](aatk::TF::ElectronStates& N, int n) -> double {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            return N(n);
        })

        .def("discrete", [](aatk::TF::ElectronStates& N) -> double {
            return N.discrete();
        })
        .def("discrete", [](aatk::TF::ElectronStates& N, double e) -> double {
            return N.discrete(e);
        })
        .def("discrete", [api](aatk::TF::ElectronStates& N, py::array_t<double> en) -> py::array {
            if (en.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions");
            }
            // get c-array representation
            auto ce   = en.data();
            auto size = en.size();
            auto cN   = N.discrete(ce, size);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          en.ndim(), 
                    (Py_intptr_t*) en.shape(), 
                    (Py_intptr_t*) en.strides(),
                    (void *)       cN, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("continuous", [](aatk::TF::ElectronStates& N) -> double {
            return N.continuous();
        })
        .def("continuous", [](aatk::TF::ElectronStates& N, double e) -> double {
            return N.continuous(e);
        })
        .def("continuous", [api](aatk::TF::ElectronStates& N, py::array_t<double> en) -> py::array {
            if (en.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions");
            }
            // get c-array representation
            auto ce   = en.data();
            auto size = en.size();
            auto cN   = N.continuous(ce, size);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          en.ndim(), 
                    (Py_intptr_t*) en.shape(), 
                    (Py_intptr_t*) en.strides(),
                    (void *)       cN, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })
        .def("setV", [](aatk::TF::ElectronStates& N, double V) {
            N.setV(V);
        })
        .def("setT", [](aatk::TF::ElectronStates& N, double T) {
            N.setT(T);
        })
        .def("setZ", [](aatk::TF::ElectronStates& N, double Z) {
            N.setZ(Z);
        })
        .def("setVTZ", [](aatk::TF::ElectronStates& N, double V, double T, double Z) {
            N.setVTZ(V,T,Z);
        })
        .def("setTolerance", [](aatk::TF::ElectronStates& N, double tol){
            N.setTolerance(tol);
        })
        .def("eBoundary", [](aatk::TF::ElectronStates& N) -> double {
            return N.eBoundary();
        })
        .def("setNmax", [](aatk::TF::ElectronStates& N, std::size_t n){
            N.setNmax(n);
        })

#ifdef ENABLE_MULTITHREADING
        .def("setThreadsLimit", [](aatk::TF::ElectronStates& N, std::size_t Nthreads) {
            N.setThreadsLimit(Nthreads);
        })
#endif
    ;

     py::class_<aatk::TF::EnergyLevel>(m, "EnergyLevel")
        .def(py::init([](){
            auto e = new aatk::TF::EnergyLevel();
            return e;
        }))

        .def("__call__", [](aatk::TF::EnergyLevel& e, int n, int l) -> double {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            if (l < 0 && l >= n) 
                throw std::runtime_error("Incorrect input: quantum number l should be between 0 and n - 1");
            return e(n, l);
        })

        .def("__getitem__", [](aatk::TF::EnergyLevel& e, int n) -> py::array {
            if (n < 1) 
                throw std::runtime_error("Incorrect input: quantum number n < 1");
            auto pyen = py::array_t<double>(n);
            auto pydata = pyen.mutable_data();
            auto en = e[n];
            for (int l = 0; l < n; ++l) pydata[l] = en[l];
            return pyen;
        })

        .def("setV", [](aatk::TF::EnergyLevel& e, double V) {
            e.setV(V);
        })
        .def("setT", [](aatk::TF::EnergyLevel& e, double T) {
            e.setT(T);
        })
        .def("setZ", [](aatk::TF::EnergyLevel& e, double Z) {
            e.setZ(Z);
        })
        .def("setVTZ", [](aatk::TF::EnergyLevel& e, double V, double T, double Z) {
            e.setVTZ(V,T,Z);
        })
        .def("setTolerance", [](aatk::TF::EnergyLevel& e, double tol){
            e.setTolerance(tol);
        })
#ifdef ENABLE_MULTITHREADING
        .def("setThreadsLimit", [](aatk::TF::EnergyLevel& e, std::size_t Nthreads) {
            e.setThreadsLimit(Nthreads);
        })
#endif
    ;

    py::class_<aatk::TF::Potential>(m, "Potential")
        .def(py::init([](){
            auto phi = new aatk::TF::Potential();
            return phi;
        }))

        .def("__call__", [](aatk::TF::Potential& phi, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return phi(x);
        })

        .def("__call__", [api](aatk::TF::Potential& phi, py::array_t<double> x) -> py::array {
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
            auto cy = phi(cx, size);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          x.ndim(), 
                    (Py_intptr_t*) x.shape(), 
                    (Py_intptr_t*) x.strides(),
                    (void *)       cy, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("dx", [](aatk::TF::Potential& phi, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return phi(x);
        })

        .def("dx", [api](aatk::TF::Potential& phi, py::array_t<double> x) -> py::array {
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
            auto cy = phi(cx, size);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          x.ndim(), 
                    (Py_intptr_t*) x.shape(), 
                    (Py_intptr_t*) x.strides(),
                    (void *)       cy, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("setV", [](aatk::TF::Potential& phi, double V) {
            phi.setV(V);
        })
        .def("setT", [](aatk::TF::Potential& phi, double T) {
            phi.setT(T);
        })
        .def("setZ", [](aatk::TF::Potential& phi, double Z) {
            phi.setZ(Z);
        })
        .def("setVTZ", [](aatk::TF::Potential& phi, double V, double T, double Z) {
            phi.setVTZ(V,T,Z);
        })
        .def("setTolerance", [](aatk::TF::Potential& phi, double tol){
            phi.setTolerance(tol);
        })
    ;

    py::class_<aatk::TF::RotatePoints>(m, "RotatePoints")
        .def(py::init([](){
            auto RP = new aatk::TF::RotatePoints();
            return RP;
        }))
        .def("inner", [](aatk::TF::RotatePoints& RP, double e, double l) -> double {
            return RP.inner(e, l);
        })
        .def("outer", [](aatk::TF::RotatePoints& RP, double e, double l) -> double {
            return RP.outer(e, l);
        })
        .def("setV", [](aatk::TF::RotatePoints& RP, double V) {
            RP.setV(V);
        })
        .def("setT", [](aatk::TF::RotatePoints& RP, double T) {
            RP.setT(T);
        })
        .def("setZ", [](aatk::TF::RotatePoints& RP, double Z) {
            RP.setZ(Z);
        })
        .def("setVTZ", [](aatk::TF::RotatePoints& RP, double V, double T, double Z) {
            RP.setVTZ(V,T,Z);
        })
        .def("setTolerance", [](aatk::TF::RotatePoints& RP, double tol) {
            RP.setTolerance(tol);
        })
    ;

}