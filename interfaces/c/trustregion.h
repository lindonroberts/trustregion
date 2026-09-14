#ifndef TRUSTREGION_H
#define TRUSTREGION_H

void f_trsunc(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    double* lambda,
    double* s
);

void f_arcunc(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    double* lambda,
    double* s
);

void f_trsapp(
    const int n, 
    const double* delta,
    const double* g,
    const double* hess,
    const double* tol,
    double* crvmin,
    double* s,
    int* info
);

void f_trsbox(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    const double* sl,
    const double* su,
    const double* tol,
    const double* xopt,
    double* crvmin,
    double* s
);

void f_trslin(
    const int n, 
    const int m, 
    const double* amat,
    const double* bvec,
    const double* xopt,
    const double delta,
    const double* g,
    const double* hess,
    const double* tol,
    double* s
);

#endif