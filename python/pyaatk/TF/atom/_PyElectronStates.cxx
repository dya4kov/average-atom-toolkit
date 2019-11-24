#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/thomas-fermi/atom/electron-states.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyElectronStates, m) {

    auto& api = py::detail::npy_api::get();

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
}