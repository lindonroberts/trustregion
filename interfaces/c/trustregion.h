/*
* The trustregion package provides a selection of algorithms for solving
* trust-region subproblems, as used in nonconvex, nonlinear optimization.
* 
* It has routines for solving trust-region subproblems of the form 
*   min_{s \in R^n} g^T * s + 0.5 * s^T * H * s
*   s.t.  ||s||_2 <= delta
* with optional bound or linear inequality constraints, or the related
* cubic regularization subproblem
*   min_{s \in R^n} g^T * s + 0.5 * s^T * H * s + (delta/3) * ||s||_2^3
* where delta >= 0 is a parameter in all cases.
* 
* For more information about trust-region and cubic regularization methods, see
* for example 
*   [CGT2000] A. R. Conn, N. I. M. Gould, P. L. Toint. Trust-Region Methods. SIAM (2000).
*   [CGT2022] C. Cartis, N. I. M. Gould, P. L. Toint. Evaluation Complexity of Algorithms for Nonconvex Optimization. SIAM (2022).
*   [NW2006] J. Nocedal, S. W. Wright. Nonlinear Optimization, 2nd ed. Springer (2006).
* 
* The implemented functions are: 
* - trsunc = globally solve the standard trust-region subproblem
* - arcunc = globally solve the cubic regularization subproblem
* - trsapp = approximately solve the standard trust-region subproblem
* - trsbox = approximately solve the bound-constrained trust-region subproblem
* - trslin = approximately solve the linear inequality-constrained trust-region subproblem
* 
* The routines trsunc and arcunc are implementations of the algorithm from
* [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
*           regularised subproblems in optimization. Mathematical Programming Computation 2:1
*           (2010), pp. 21-57. 
* 
* The routines trsapp, trsbox and trslin are implementations of the CG-Steihaug method
* (with constraints handled via an active set approach) copied with minimal modifications
* from the PRIMA software of Z. Zhang: https://github.com/libprima/prima/tree/main
* They are improved versions of the subproblem solvers from M. J. D. Powell's 
* NEWUOA, BOBYQA and LINCOA codes.
* 
* PRIMA is made available under the following BSD-3 clause licence:
* (see https://github.com/libprima/prima/blob/main/LICENCE.txt)
* 
* PRIMA BSD 3-Clause License
* 
* Copyright (c) 2020--2026, Zaikun ZHANG ( https://www.zhangzk.net )
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* 
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 
* 3. Neither the name of the copyright holder nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TRUSTREGION_H
#define TRUSTREGION_H

#ifdef __cplusplus
extern "C" {
#endif

/*
* Globally solve the (unconstrained) trust-region subproblem
*   min_{s \in R^n} g^T * s + 0.5 * s^T * H * s
*   s.t.  ||s||_2 <= delta
*
* This algorithm is largely suited for small/medium scale problems
* (i.e. n not too large, e.g. n <= 500)
* 
* Inputs are:
* - n = dimension of problem
* - delta = non-negative trust-region radius
* - g = vector of length n
* - hess = array of length n*n storing symmetric H in row/column-major ordering
*          (both orderings are equivalent for symmetric matrices)
* - lambda = Lagrange multiplier at solution, >= 0
* - s = global minimizer of length n
* 
* This function implements the algorithm from
* [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
*           regularised subproblems in optimization. Mathematical Programming Computation 2:1
*           (2010), pp. 21-57. 
*/
void f_trsunc(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    double* lambda,
    double* s
);

/*
* Globally solve the (unconstrained) cubic regularization subproblem
*   min_{s \in R^n} g^T * s + 0.5 * s^T * H * s + (delta/3) * ||s||_2^3
* 
* This algorithm is largely suited for small/medium scale problems
* (i.e. n not too large, e.g. n <= 500)
* 
* Inputs are:
* - n = dimension of problem
* - delta = non-negative regularization parameter
* - g = vector of length n
* - hess = array of length n*n storing symmetric H in row/column-major ordering
*          (both orderings are equivalent for symmetric matrices)
* - s = global minimizer of length n
* 
* This function implements the algorithm from
* [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
*           regularised subproblems in optimization. Mathematical Programming Computation 2:1
*           (2010), pp. 21-57. 
*/
void f_arcunc(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    double* s
);

/*
* Approximately solve the (unconstrained) trust-region subproblem
*   min_{s \in R^n} g^T * s + 0.5 * s^T * H * s
*   s.t.  ||s||_2 <= delta
* 
* Inputs are:
* - n = dimension of problem
* - delta = non-negative trust-region radius
* - g = vector of length n
* - hess = array of length n*n storing symmetric H in row/column-major ordering
*          (both orderings are equivalent for symmetric matrices)
* - tol = strictly positive termination tolerance
* - crvmin = estimate of minimum eigenvalue of H (based on observed iterates)
* - s = approximate global minimizer of length n
* - info = output information flag
*          info = 0 --> approximate solution found
*          info = 1 --> last iteration gave insufficient objective reduction
*          info = 2 --> max iterations reached
*          info = -1 --> rounding errors too large
* 
* This function implements a modified version of the CG-Steihaug method
* from https://github.com/libprima/prima/blob/main/fortran/newuoa/trustregion.f90
* (an improvement of the original code from M. J. D. Powell's NEWUOA software [Pow2006])
* 
* [Pow2006] M. J. D. Powell, The NEWUOA Software for Unconstrained Optimization without Derivatives,
*           in: Large-Scale Nonlinear Optimization, P. Pardalos, G. Di Pillo, M. Roma, eds,
*           Springer, 2006, pp. 225-297.
*/
void f_trsapp(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    const double tol,
    double* crvmin,
    double* s,
    int* info
);

/*
* Approximately solve the box-constrained trust-region subproblem
*   min_{s \in R^n} g^T * s + 0.5 * s^T * H * s
*   s.t.  ||s||_2 <= delta
*         sl <= xopt + s <= su
* 
* The code assumes that the point s=0 is feasible (i.e. sl <= xopt <= su)
* 
* Inputs are:
* - n = dimension of problem
* - delta = non-negative trust-region radius
* - g = vector of length n
* - hess = array of length n*n storing symmetric H in row/column-major ordering
*          (both orderings are equivalent for symmetric matrices)
* - sl = lower bounds, array of length n
* - su = upper bounds, array of length n
* - xopt = base point for bound constraints, array of length n
* - tol = strictly positive termination tolerance
* - crvmin = estimate of minimum eigenvalue of H (based on observed iterates)
* - s = approximate global minimizer of length n
* 
* This function implements a modified version of the active set/CG-Steihaug method
* from https://github.com/libprima/prima/blob/main/fortran/bobyqa/trustregion.f90
* (an improvement of the original code from M. J. D. Powell's BOBYQA software [Pow2009])
* 
* [Pow2009] M. J. D. Powell, The BOBYQA Algorithm for Bound Constrained Optimization without Derivatives,
*           Technical report DAMTP 2009/NA06, University of Cambridge, 2006.
*/
void f_trsbox(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    const double* sl,
    const double* su,
    const double* xopt,
    const double tol,
    double* crvmin,
    double* s
);

/*
* Approximately solve the linearly-constrained trust-region subproblem
*   min_{s \in R^n} g^T * s + 0.5 * s^T * H * s
*   s.t.  ||s||_2 <= delta
*         A * (xopt + s) <= bvec
* 
* The code assumes that the point s=0 is feasible (i.e. A * xopt <= bvec)
* 
* Inputs are:
* - n = dimension of problem
* - delta = non-negative trust-region radius
* - g = vector of length n
* - hess = array of length n*n storing symmetric H in row/column-major ordering
*          (both orderings are equivalent for symmetric matrices)
* - m = number of linear inequality constraints
* - amat = constraint LHS matrix, array of length m*n in C-style memory layout (row-major ordering)
*          i.e. A(0, :) = [amat[0], ..., amat[n-1]], and A(1, :) = [amat[n], ..., amat[2*n-1]], etc.
*          or A(i, j) = amat[i * n + j] for row index i=0, ..., m-1 and column index j=0, ..., n-1
* - bvec = constraint RHS vector, array of length m
* - xopt = base point for constraint, array of length n
* - tol = strictly positive termination tolerance
* - s = approximate global minimizer of length n
* 
* This function implements a modified version of the active set/CG-Steihaug method
* from https://github.com/libprima/prima/blob/main/fortran/lincoa/trustregion.f90
* (an improvement of the original code from M. J. D. Powell's LINCOA software [Pow2015])
* 
* [Pow2015] M. J. D. Powell, On Fast Trust Region Methods for Quadratic Models with Linear Constraints,
*           Mathematical Programming Computation 7:3 (2015), pp. 237-267.
*/
void f_trslin(
    const int n, 
    const double delta,
    const double* g,
    const double* hess,
    const int m, 
    const double* amat,
    const double* bvec,
    const double* xopt,
    const double tol,
    double* s
);

#ifdef __cplusplus
}
#endif

#endif