#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "trustregion.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

typedef py::array_t<double> NumpyNArray;

__attribute__((visibility("default"))) 
py::tuple py_trsunc(double delta, NumpyNArray g, NumpyNArray H)
{
    py::buffer_info buf_g = g.request();
    py::buffer_info buf_H = H.request();
 
    // Check dimensions of all inputs
    if (buf_g.ndim != 1) throw std::runtime_error("g must be a vector");
    if (buf_H.ndim != 2) throw std::runtime_error("H must be a matrix");

    int n;
    n = (int)buf_g.size;
    
    if (buf_H.shape[0] != n) throw std::runtime_error("H has incompatible number of rows with g");
    if (buf_H.shape[1] != n) throw std::runtime_error("H has incompatible number of columns with g");
    
    double* ptr_g = static_cast<double*>(buf_g.ptr);
    double* ptr_H = static_cast<double*>(buf_H.ptr);
    
    // Define output info
    double lambda;
    NumpyNArray s = NumpyNArray(n);
    py::buffer_info buf_s = s.request();
    double* ptr_s = static_cast<double*>(buf_s.ptr);

    f_trsunc(n, delta, ptr_g, ptr_H, &lambda, ptr_s);
    
    // Returns tuple (s, lambda)
    py::tuple output = py::make_tuple(s, lambda);
    return output;
}

__attribute__((visibility("default"))) 
py::tuple py_arcunc(double delta, NumpyNArray g, NumpyNArray H)
{
    py::buffer_info buf_g = g.request();
    py::buffer_info buf_H = H.request();
 
    // Check dimensions of all inputs
    if (buf_g.ndim != 1) throw std::runtime_error("g must be a vector");
    if (buf_H.ndim != 2) throw std::runtime_error("H must be a matrix");

    int n;
    n = (int)buf_g.size;
    
    if (buf_H.shape[0] != n) throw std::runtime_error("H has incompatible number of rows with g");
    if (buf_H.shape[1] != n) throw std::runtime_error("H has incompatible number of columns with g");
    
    double* ptr_g = static_cast<double*>(buf_g.ptr);
    double* ptr_H = static_cast<double*>(buf_H.ptr);
    
    // Define output info
    NumpyNArray s = NumpyNArray(n);
    py::buffer_info buf_s = s.request();
    double* ptr_s = static_cast<double*>(buf_s.ptr);

    f_arcunc(n, delta, ptr_g, ptr_H, ptr_s);
    
    // Returns tuple (s,)
    py::tuple output = py::make_tuple(s);
    return output;
}

__attribute__((visibility("default"))) 
py::tuple py_trsapp(double delta, NumpyNArray g, NumpyNArray H, double tol)
{
    py::buffer_info buf_g = g.request();
    py::buffer_info buf_H = H.request();
 
    // Check dimensions of all inputs
    if (buf_g.ndim != 1) throw std::runtime_error("g must be a vector");
    if (buf_H.ndim != 2) throw std::runtime_error("H must be a matrix");

    int n;
    n = (int)buf_g.size;
    
    if (buf_H.shape[0] != n) throw std::runtime_error("H has incompatible number of rows with g");
    if (buf_H.shape[1] != n) throw std::runtime_error("H has incompatible number of columns with g");
    
    double* ptr_g = static_cast<double*>(buf_g.ptr);
    double* ptr_H = static_cast<double*>(buf_H.ptr);
    
    // Define output info
    double crvmin;
    NumpyNArray s = NumpyNArray(n);
    py::buffer_info buf_s = s.request();
    double* ptr_s = static_cast<double*>(buf_s.ptr);
    int info;

    f_trsapp(n, delta, ptr_g, ptr_H, tol, &crvmin, ptr_s, &info);
    
    // Returns tuple (s,)
    py::tuple output = py::make_tuple(s, crvmin, info);
    return output;
}

PYBIND11_MODULE(_trustregion, m) {
    m.doc() = R"pbdoc(
        trustregion simple interface
        ----------------------

        .. currentmodule:: _trustregion

        .. autosummary::
           :toctree: _generate

           py_trsunc
           py_arcunc
           py_trsapp
    )pbdoc";
    
    m.def("py_trsunc", &py_trsunc, R"pbdoc(
        Call trsunc solver
    )pbdoc");

    m.def("py_arcunc", &py_arcunc, R"pbdoc(
        Call arcunc solver
    )pbdoc");

    m.def("py_trsapp", &py_trsapp, R"pbdoc(
        Call trsapp solver
    )pbdoc");


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
