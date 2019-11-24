#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyChemicalPotential, m) {

    auto& api = py::detail::npy_api::get();

    py::class_<aatk::TF::ChemicalPotential>(m, "ChemicalPotential")
        .def(py::init([](){
            auto M = new aatk::TF::ChemicalPotential();
            return M;
        }))

        .def("__call__", [](aatk::TF::ChemicalPotential& M, double V, double T) -> double {
            return M(V,T);
        })

        .def("__call__", [api](aatk::TF::ChemicalPotential& M, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Mdata = M(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("__call__", [api](aatk::TF::ChemicalPotential& M, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Mdata = M(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("__call__", [api](aatk::TF::ChemicalPotential& M, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Mdata = M(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DV", [](aatk::TF::ChemicalPotential& M, double V, double T) -> double {
            return M.DV(V,T);
        })

        .def("DV", [api](aatk::TF::ChemicalPotential& M, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Mdata = M.DV(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DV", [api](aatk::TF::ChemicalPotential& M, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Mdata = M.DV(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DV", [api](aatk::TF::ChemicalPotential& M, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Mdata = M.DV(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DT", [](aatk::TF::ChemicalPotential& M, double V, double T) -> double {
            return M.DT(V,T);
        })

        .def("DT", [api](aatk::TF::ChemicalPotential& M, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Mdata = M.DT(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DT", [api](aatk::TF::ChemicalPotential& M, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Mdata = M.DT(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DT", [api](aatk::TF::ChemicalPotential& M, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Mdata = M.DT(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Mdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("setZ", [](aatk::TF::ChemicalPotential& M, double Z) {
            M.setZ(Z);
        })
        .def("setTolerance", [](aatk::TF::ChemicalPotential& M, double tol){
            M.setTolerance(tol);
        })
#ifdef ENABLE_MULTITHREADING
        .def("setThreadsLimit", [](aatk::TF::ChemicalPotential& M, std::size_t Nthreads) {
            M.setThreadsLimit(Nthreads);
        })
#endif
    ;
}