#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/thomas-fermi/atom/shell/potential.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyPotentialShell, m) {

    auto& api = py::detail::npy_api::get();

    py::class_<aatk::TF::shell::Potential>(m, "Potential")
        .def(py::init([](){
            auto phi = new aatk::TF::shell::Potential();
            return phi;
        }))

        .def("__call__", [](aatk::TF::shell::Potential& phi, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return phi(x);
        })
        .def("__call__", [api](aatk::TF::shell::Potential& phi, py::array_t<double> x) -> py::array {
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

        .def("dx", [](aatk::TF::shell::Potential& phi, double x) -> double {
            bool correct_input = (x >= 0.0 && x <= 1.0);
            if (!correct_input) {
                throw std::runtime_error("Incorrect input: x should be between 0 and 1");
            }
            return phi(x);
        })

        .def("dx", [api](aatk::TF::shell::Potential& phi, py::array_t<double> x) -> py::array {
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
        .def("setNmax", [](aatk::TF::shell::Potential& phi, std::size_t nmax) {
            phi.setNmax(nmax);
        })
        .def("setV", [](aatk::TF::shell::Potential& phi, double V) {
            phi.setV(V);
        })
        .def("setT", [](aatk::TF::shell::Potential& phi, double T) {
            phi.setT(T);
        })
        .def("setZ", [](aatk::TF::shell::Potential& phi, double Z) {
            phi.setZ(Z);
        })
        .def("setVTZ", [](aatk::TF::shell::Potential& phi, double V, double T, double Z) {
            phi.setVTZ(V,T,Z);
        })
        .def("setTolerance", [](aatk::TF::shell::Potential& phi, double tol){
            phi.setTolerance(tol);
        })
    ;
}