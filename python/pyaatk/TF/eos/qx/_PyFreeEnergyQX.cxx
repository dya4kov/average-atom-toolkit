#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/thomas-fermi/eos/qx/free-energy.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyFreeEnergyQX, m) {

    auto& api = py::detail::npy_api::get();

    py::class_<aatk::TF::qx::FreeEnergy>(m, "FreeEnergy")
        .def(py::init([](){
            auto F = new aatk::TF::qx::FreeEnergy();
            return F;
        }))

        .def("__call__", [](aatk::TF::qx::FreeEnergy& F, double V, double T) -> double {
            return F(V,T);
        })

        .def("__call__", [api](aatk::TF::qx::FreeEnergy& F, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Fdata = F(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("__call__", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Fdata = F(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("__call__", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Fdata = F(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DV", [](aatk::TF::qx::FreeEnergy& F, double V, double T) -> double {
            return F.DV(V,T);
        })

        .def("DV", [api](aatk::TF::qx::FreeEnergy& F, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Fdata = F.DV(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DV", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Fdata = F.DV(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DV", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Fdata = F.DV(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DT", [](aatk::TF::qx::FreeEnergy& F, double V, double T) -> double {
            return F.DT(V,T);
        })

        .def("DT", [api](aatk::TF::qx::FreeEnergy& F, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Fdata = F.DT(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DT", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Fdata = F.DT(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DT", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Fdata = F.DT(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("D2V", [](aatk::TF::qx::FreeEnergy& F, double V, double T) -> double {
            return F.D2V(V,T);
        })

        .def("D2V", [api](aatk::TF::qx::FreeEnergy& F, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Fdata = F.D2V(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("D2V", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Fdata = F.D2V(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("D2V", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Fdata = F.D2V(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DVT", [](aatk::TF::qx::FreeEnergy& F, double V, double T) -> double {
            return F.DVT(V,T);
        })

        .def("DVT", [api](aatk::TF::qx::FreeEnergy& F, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Fdata = F.DVT(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DVT", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Fdata = F.DVT(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("DVT", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Fdata = F.DVT(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("D2T", [](aatk::TF::qx::FreeEnergy& F, double V, double T) -> double {
            return F.D2T(V,T);
        })

        .def("D2T", [api](aatk::TF::qx::FreeEnergy& F, double V, py::array_t<double> T) -> py::array {
            if (T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Temperature should be 1D array");
            }
            // get c-array representation
            const double* Vdata = &V;
            const double* Tdata = T.data();
            std::size_t tsize = T.size();
            std::size_t vsize = 1;
            auto Fdata = F.D2T(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          T.ndim(), 
                    (Py_intptr_t*) T.shape(), 
                    (Py_intptr_t*) T.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("D2T", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, double T) -> py::array {
            if (V.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. Volume should be 1D array");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = &T;
            std::size_t vsize = V.size();
            std::size_t tsize = 1;
            auto Fdata = F.D2T(Vdata, Tdata, vsize, tsize);

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          V.ndim(), 
                    (Py_intptr_t*) V.shape(), 
                    (Py_intptr_t*) V.strides(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("D2T", [api](aatk::TF::qx::FreeEnergy& F, py::array_t<double> V, py::array_t<double> T) -> py::array {
            if (V.ndim() != 1 || T.ndim() != 1) {
                throw std::runtime_error("Incorrect number of dimensions. V and T should be 1D arrays");
            }
            // get c-array representation
            const double* Vdata = V.data();
            const double* Tdata = T.data();
            std::size_t vsize = V.size();
            std::size_t tsize = T.size();
            auto Fdata = F.D2T(Vdata, Tdata, vsize, tsize);
            std::vector<std::size_t> shape({vsize, tsize});
            std::vector<std::size_t> strides({tsize*sizeof(double), sizeof(double)});

            return py::reinterpret_steal<py::array>(
                api.PyArray_NewFromDescr_(
                    api.PyArray_Type_, 
                    py::dtype::of<double>().release().ptr(),
                    (int)          2, 
                    (Py_intptr_t*) shape.data(), 
                    (Py_intptr_t*) strides.data(),
                    (void *)       Fdata, 
                    py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                    py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                    py::detail::npy_api::NPY_ARRAY_OWNDATA_ |
                    py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                    nullptr
                )
            );
        })

        .def("setZ", [](aatk::TF::qx::FreeEnergy& F, double Z) {
            F.setZ(Z);
        })
        .def("setTolerance", [](aatk::TF::qx::FreeEnergy& F, double tol){
            F.setTolerance(tol);
        })
#ifdef ENABLE_MULTITHREADING
        .def("setThreadsLimit", [](aatk::TF::qx::FreeEnergy& F, std::size_t Nthreads) {
            F.setThreadsLimit(Nthreads);
        })
#endif
    ;
}