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
NumpyNArray py_arcunc(double delta, NumpyNArray g, NumpyNArray H)
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
    
    return s;
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
    
    // Returns tuple (s, crvmin, info)
    py::tuple output = py::make_tuple(s, crvmin, info);
    return output;
}

__attribute__((visibility("default"))) 
py::tuple py_trsbox(double delta, NumpyNArray g, NumpyNArray H, NumpyNArray sl, NumpyNArray su, NumpyNArray xopt, double tol)
{
    py::buffer_info buf_g = g.request();
    py::buffer_info buf_H = H.request();
    py::buffer_info buf_sl = sl.request();
    py::buffer_info buf_su = su.request();
    py::buffer_info buf_xopt = xopt.request();
 
    // Check dimensions of all inputs
    if (buf_g.ndim != 1) throw std::runtime_error("g must be a vector");
    if (buf_H.ndim != 2) throw std::runtime_error("H must be a matrix");
    if (buf_sl.ndim != 1) throw std::runtime_error("sl must be a vector");
    if (buf_su.ndim != 1) throw std::runtime_error("su must be a vector");
    if (buf_xopt.ndim != 1) throw std::runtime_error("xopt must be a vector");

    int n;
    n = (int)buf_g.size;
    
    if (buf_H.shape[0] != n) throw std::runtime_error("H has incompatible number of rows with g");
    if (buf_H.shape[1] != n) throw std::runtime_error("H has incompatible number of columns with g");
    if (buf_sl.size != n) throw std::runtime_error("sl has incompatible dimensions with g");
    if (buf_su.size != n) throw std::runtime_error("su has incompatible dimensions with g");
    if (buf_xopt.size != n) throw std::runtime_error("xopt has incompatible dimensions with g");
    
    double* ptr_g = static_cast<double*>(buf_g.ptr);
    double* ptr_H = static_cast<double*>(buf_H.ptr);
    double* ptr_sl = static_cast<double*>(buf_sl.ptr);
    double* ptr_su = static_cast<double*>(buf_su.ptr);
    double* ptr_xopt = static_cast<double*>(buf_xopt.ptr);
    
    // Define output info
    double crvmin;
    NumpyNArray s = NumpyNArray(n);
    py::buffer_info buf_s = s.request();
    double* ptr_s = static_cast<double*>(buf_s.ptr);

    f_trsbox(n, delta, ptr_g, ptr_H, ptr_sl, ptr_su, ptr_xopt, tol, &crvmin, ptr_s);
    
    // Returns tuple (s, crvmin)
    py::tuple output = py::make_tuple(s, crvmin);
    return output;
}

__attribute__((visibility("default"))) 
NumpyNArray py_trslin(double delta, NumpyNArray g, NumpyNArray H, NumpyNArray A, NumpyNArray b, NumpyNArray xopt, double tol)
{
    py::buffer_info buf_g = g.request();
    py::buffer_info buf_H = H.request();
    py::buffer_info buf_A = A.request();
    py::buffer_info buf_b = b.request();
    py::buffer_info buf_xopt = xopt.request();
 
    // Check dimensions of all inputs
    if (buf_g.ndim != 1) throw std::runtime_error("g must be a vector");
    if (buf_H.ndim != 2) throw std::runtime_error("H must be a matrix");
    if (buf_A.ndim != 2) throw std::runtime_error("A must be a matrix");
    if (buf_b.ndim != 1) throw std::runtime_error("b must be a vector");
    if (buf_xopt.ndim != 1) throw std::runtime_error("xopt must be a vector");

    int n;
    int m;
    n = (int)buf_g.size;
    m = (int)buf_A.shape[0];

    // If no linear constraints, then A = [[]] has size (1,0)
    if (buf_A.shape[0] * buf_A.shape[1] == 0) m = 0;
    
    if (buf_H.shape[0] != n) throw std::runtime_error("H has incompatible number of rows with g");
    if (buf_H.shape[1] != n) throw std::runtime_error("H has incompatible number of columns with g");
    if ((m > 0) && (buf_A.shape[1] != n)) throw std::runtime_error("A has incompatible number of columns with g");
    if (buf_b.size != m) throw std::runtime_error("b has incompatible dimensions with A");
    if (buf_xopt.size != n) throw std::runtime_error("xopt has incompatible dimensions with g");
    
    double* ptr_g = static_cast<double*>(buf_g.ptr);
    double* ptr_H = static_cast<double*>(buf_H.ptr);
    double* ptr_A = static_cast<double*>(buf_A.ptr);
    double* ptr_b = static_cast<double*>(buf_b.ptr);
    double* ptr_xopt = static_cast<double*>(buf_xopt.ptr);
    
    // Define output info
    double crvmin;
    NumpyNArray s = NumpyNArray(n);
    py::buffer_info buf_s = s.request();
    double* ptr_s = static_cast<double*>(buf_s.ptr);

    f_trslin(n, delta, ptr_g, ptr_H, m, ptr_A, ptr_b, ptr_xopt, tol, ptr_s);
    
    return s;
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

    m.def("py_trsbox", &py_trsbox, R"pbdoc(
        Call trsbox solver
    )pbdoc");

    m.def("py_trslin", &py_trslin, R"pbdoc(
        Call trslin solver
    )pbdoc");


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
