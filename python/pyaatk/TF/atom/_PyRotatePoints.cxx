#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/thomas-fermi/atom/rotate-points.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyRotatePoints, m) {

    auto& api = py::detail::npy_api::get();

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