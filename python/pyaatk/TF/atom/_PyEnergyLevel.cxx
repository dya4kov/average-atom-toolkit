#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <average-atom-toolkit/thomas-fermi/atom/energy-level.h>

namespace py = pybind11;

PYBIND11_MODULE(_PyEnergyLevel, m) {

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
}