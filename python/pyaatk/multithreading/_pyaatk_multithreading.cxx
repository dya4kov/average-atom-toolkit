#include <average-atom-toolkit/multithreading/thread-pool.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

using ::aatk::multithreading::ThreadPool;

PYBIND11_MODULE(_pyaatk_multithreading, m) {

py::class_<ThreadPool>(m, "ThreadPool")
    .def(py::init([](std::size_t nthreads) {
        auto pool = new ThreadPool(nthreads);
        return pool;
    }),
        py::arg("count") = 0
    )
    .def_property_readonly("size", [](ThreadPool& pool) -> std::size_t {
        return pool.size();
    })
    .def("resize", [](ThreadPool& pool, std::size_t nthreads) -> void {
        pool.resize(nthreads);
    })
    .def("start", [](ThreadPool& pool) -> void {
        pool.start();
    })
    .def("stop", [](ThreadPool& pool) -> void {
        pool.stop();
    });
}
