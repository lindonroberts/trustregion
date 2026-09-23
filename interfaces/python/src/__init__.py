"""
The trustregion package provides a selection of algorithms for solving
trust-region subproblems, as used in nonconvex, nonlinear optimization.

It has routines for solving trust-region subproblems of the form 
  min_{s in R^n} g^T * s + 0.5 * s^T * H * s
  s.t.  ||s||_2 <= delta
with optional bound or linear inequality constraints, or the related
cubic regularization subproblem
  min_{s in R^n} g^T * s + 0.5 * s^T * H * s + (delta/3) * ||s||_2^3
where delta >= 0 is a parameter in all cases.

For more information about trust-region and cubic regularization methods, see
for example 
  [CGT2000] A. R. Conn, N. I. M. Gould, P. L. Toint. Trust-Region Methods. SIAM (2000).
  [CGT2022] C. Cartis, N. I. M. Gould, P. L. Toint. Evaluation Complexity of Algorithms for Nonconvex Optimization. SIAM (2022).
  [NW2006] J. Nocedal, S. W. Wright. Nonlinear Optimization, 2nd ed. Springer (2006).

The implemented functions are: 
- trsunc = globally solve the standard trust-region subproblem
- arcunc = globally solve the cubic regularization subproblem
- trsapp = approximately solve the standard trust-region subproblem
- trsbox = approximately solve the bound-constrained trust-region subproblem
- trslin = approximately solve the linear inequality-constrained trust-region subproblem

The routines trsunc and arcunc are implementations of the algorithm from
[GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
          regularised subproblems in optimization. Mathematical Programming Computation 2:1
          (2010), pp. 21-57. 

The routines trsapp, trsbox and trslin are implementations of the CG-Steihaug method
(with constraints handled via an active set approach) copied with minimal modifications
from the PRIMA software of Z. Zhang: https://github.com/libprima/prima/tree/main
They are improved versions of the subproblem solvers from M. J. D. Powell's 
NEWUOA, BOBYQA and LINCOA codes.

PRIMA is made available under the following BSD-3 clause licence:
(see https://github.com/libprima/prima/blob/main/LICENCE.txt)

PRIMA BSD 3-Clause License

Copyright (c) 2020--2026, Zaikun ZHANG ( https://www.zhangzk.net )

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

__version__ = '0.1.0'

import numpy as np
from ._trustregion import py_trsunc, py_arcunc, py_trsapp, py_trsbox, py_trslin

__all__ = ['trsunc', 'arcunc', 'trsapp', 'trsbox', 'trslin']


def force_1d_array(x, name):
    """
    Force x to be a 1-D numpy array of floating-point type.
    
    Input 'name' is a string to be used when printing an error message
    """
    x = np.atleast_1d(np.asarray(x))

    if x.ndim != 1:
        raise ValueError("'%s' must only have one dimension." % str(name))

    if x.dtype.kind in np.typecodes["AllInteger"]:
        x = np.asarray(x, dtype=float)
    
    return x


def force_2d_array(A, name):
    """
    Force A to be a 2-D numpy array of floating-point type.
    
    Input 'name' is a string to be used when printing an error message
    """
    A = np.atleast_2d(np.asarray(A))

    if A.ndim != 2:
        raise ValueError("'%s' must have exactly two dimensions." % str(name))

    if A.dtype.kind in np.typecodes["AllInteger"]:
        A = np.asarray(A, dtype=float)
    
    return A


def validate_inputs(delta, g, H, sl=None, su=None, xopt=None, A=None, b=None, thresh=1e-14):
    if delta < 0.0:
        raise RuntimeError("delta must be non-negative")

    g = force_1d_array(g, "g")
    H = force_2d_array(H, "H")

    # Check dimensions work out
    n = len(g)
    if H.shape[0] != n:
        raise RuntimeError("H has number of rows incompatible with g")
    if H.shape[1] != n:
        raise RuntimeError("H has number of columns incompatible with g")

    # Check data
    if not np.all(np.isfinite(g)):
        raise RuntimeError("g cannot have nan/inf values")
    if not np.all(np.isfinite(H)):
        raise RuntimeError("H cannot have nan/inf values")
    if not np.allclose(H, H.T):
        raise RuntimeError("H must be symmetric")

    # Check presence of constraints make sense
    if (sl is not None and su is None) or (sl is None and su is not None):
        raise RuntimeError("Must specify either none or both of sl and su")
    if (A is not None and b is None) or (A is None and b is not None):
        raise RuntimeError("Must specify either none or both of A and b")
    if (sl is not None and A is not None):
        raise RuntimeError("Must specify either bounds or linear constraints, not both")
    if (sl is not None or A is not None) and xopt is None:
        raise RuntimeError("Must specify xopt if either bounds or linear constraints provided")

    # Identify constraint type
    if sl is not None:
        cons_type = 'bounds'
    elif A is not None:
        cons_type = 'linear'
    else:
        cons_type = 'none'

    if cons_type != 'none':
        xopt = force_1d_array(xopt, "xopt")
        if len(xopt) != n:
            raise RuntimeError("xopt has incompatible shape with g")
        if not np.all(np.isfinite(xopt)):
            raise RuntimeError("xopt cannot have nan/inf values")
    
    if cons_type == 'bounds':
        sl = force_1d_array(sl, "sl")
        sl[np.isnan(sl)] = -np.inf
        if len(sl) != n:
            raise RuntimeError("sl has incompatible shape with g")
        
        su = force_1d_array(su, "su")
        su[np.isnan(su)] = np.inf
        if len(su) != n:
            raise RuntimeError("su has incompatible shape with g")

        # Check xopt is feasible up to rounding errors
        if not np.all(xopt >= sl - thresh):
            raise RuntimeError("Must have xopt >= sl")
        if not np.all(xopt <= su + thresh):
            raise RuntimeError("Must have xopt <= su")

    elif cons_type == 'linear':
        A = force_2d_array(A, "A")
        m = A.shape[0]
        if A.shape[1] != n:
            raise RuntimeError("A has number of columns incompatible with g")
        
        b = force_1d_array(b, "b")
        b[np.isnan(b)] = np.inf

        if len(b) != m:
            raise RuntimeError("b has incompatible shape with A")

        # Check xopt is feasible, up to rounding errors
        if not np.all(A @ xopt <= b + thresh):
            raise RuntimeError("Must have A @ xopt <= b")

    return delta, g, H, sl, su, xopt, A, b, cons_type


def trsunc(delta, g, H):
    """
    Globally solve the (unconstrained) trust-region subproblem
        min_{s in R^n} g^T * s + 0.5 * s^T * H * s
        s.t.  ||s||_2 <= delta
    
    This algorithm is largely suited for small/medium scale problems
    (i.e. n not too large, e.g. n <= 500)

    Call with:
        s, lda = trustregion.trsunc(delta, g, H)
    
    Inputs are:
    - delta = non-negative trust-region radius
    - g = np.ndarray of shape (n,)
    - H = symmetric np.ndarray of shape (n,n)

    Outputs are:
    - s = global minimizer, np.ndarray of shape (n,)
    - lda = Lagrange multiplier at solution, >= 0
    
    This function implements the algorithm from
    [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
                regularised subproblems in optimization. Mathematical Programming Computation 2:1
                (2010), pp. 21-57. 
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H)
    s, lda = py_trsunc(delta, g, H)
    return s, lda


def arcunc(delta, g, H):
    """
    Globally solve the (unconstrained) cubic regularization subproblem
        min_{s in R^n} g^T * s + 0.5 * s^T * H * s + (delta/3) * ||s||_2^3
    
    This algorithm is largely suited for small/medium scale problems
    (i.e. n not too large, e.g. n <= 500)

    Call with:
        s = trustregion.arcunc(delta, g, H)
    
    Inputs are:
    - delta = non-negative trust-region radius
    - g = np.ndarray of shape (n,)
    - H = symmetric np.ndarray of shape (n,n)

    Output is:
    - s = global minimizer, np.ndarray of shape (n,)
    
    This function implements the algorithm from
    [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
              regularised subproblems in optimization. Mathematical Programming Computation 2:1
              (2010), pp. 21-57. 
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H)
    s = py_arcunc(delta, g, H)
    return s


def trsapp(delta, g, H, tol=1e-2):
    """
    Approximately solve the (unconstrained) trust-region subproblem
        min_{s in R^n} g^T * s + 0.5 * s^T * H * s
        s.t.  ||s||_2 <= delta

    Call with:
        s, crvmin, info = trustregion.trsapp(delta, g, H, tol=1e-2)
    
    Inputs are:
    - delta = non-negative trust-region radius
    - g = np.ndarray of shape (n,)
    - H = symmetric np.ndarray of shape (n,n)
    - tol = strictly positive termination tolerance

    Outputs are:
    - s = approximate minimizer, np.ndarray of shape (n,)
    - crvmin = estimate of minimum eigenvalue of H (based on observed iterates)
    - info = output information flag
        info = 0 --> approximate solution found
        info = 1 --> last iteration gave insufficient objective reduction
        info = 2 --> max iterations reached
        info = -1 --> rounding errors too large
    
    This function implements a modified version of the CG-Steihaug method
    from https://github.com/libprima/prima/blob/main/fortran/newuoa/trustregion.f90
    (an improvement of the original code from M. J. D. Powell's NEWUOA software [Pow2006])
    
    [Pow2006] M. J. D. Powell, The NEWUOA Software for Unconstrained Optimization without Derivatives,
              in: Large-Scale Nonlinear Optimization, P. Pardalos, G. Di Pillo, M. Roma, eds,
              Springer, 2006, pp. 225-297.
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H)
    if tol <= 0.0:
        raise RuntimeError("tol must be strictly positive")

    s, crvmin, info = py_trsapp(delta, g, H, tol)
    return s, crvmin, info


def trsbox(delta, g, H, sl, su, xopt, tol=1e-2):
    """
    Approximately solve the box-constrained trust-region subproblem
        min_{s in R^n} g^T * s + 0.5 * s^T * H * s
        s.t.  ||s||_2 <= delta
              sl <= xopt + s <= su
    
    The code assumes that the point s=0 is feasible (i.e. sl <= xopt <= su)

    Call with
        s, crvmin = trustregion.trsbox(delta, g, H, sl, su, xopt, tol=1e-2)
    
    Inputs are:
    - delta = non-negative trust-region radius
    - g = np.ndarray of shape (n,)
    - H = symmetric np.ndarray of shape (n,n)
    - sl = lower bounds, np.ndarray of shape (n,)
    - su = upper bounds, np.ndarray of shape (n,)
    - xopt = base point for bound constraints, np.ndarray of shape (n,)
    - tol = strictly positive termination tolerance

    Outputs are:
    - s = approximate minimizer, np.ndarray of shape (n,)
    - crvmin = estimate of minimum eigenvalue of H (based on observed iterates)
    
    This function implements a modified version of the active set/CG-Steihaug method
    from https://github.com/libprima/prima/blob/main/fortran/bobyqa/trustregion.f90
    (an improvement of the original code from M. J. D. Powell's BOBYQA software [Pow2009])
    
    [Pow2009] M. J. D. Powell, The BOBYQA Algorithm for Bound Constrained Optimization without Derivatives,
              Technical report DAMTP 2009/NA06, University of Cambridge, 2006.
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H, sl=sl, su=su, xopt=xopt)
    if tol <= 0.0:
        raise RuntimeError("tol must be strictly positive")

    s, crvmin = py_trsbox(delta, g, H, sl, su, xopt, tol)
    return s, crvmin


def trslin(delta, g, H, A, b, xopt, tol=1e-2):
    """
    Approximately solve the linearly-constrained trust-region subproblem
        min_{s in R^n} g^T * s + 0.5 * s^T * H * s
        s.t.  ||s||_2 <= delta
              A @ (xopt + s) <= b
    
    The code assumes that the point s=0 is feasible (i.e. A @ xopt <= b)

    Call with:
        s = trustregion.trslin(delta, g, H, A, b, xopt, tol=1e-2)
    
    Inputs are:
    - delta = non-negative trust-region radius
    - g = np.ndarray of shape (n,)
    - H = symmetric np.ndarray of shape (n,n)
    - A = constraint matrix, np.ndarray of shape (m,n)
    - b = constraint RHS vector, np.ndarray of shape (m,)
    - xopt = base point for bound constraints, np.ndarray of shape (n,)
    - tol = strictly positive termination tolerance

    Output is:
    - s = approximate minimizer, np.ndarray of shape (n,)
    
    This function implements a modified version of the active set/CG-Steihaug method
    from https://github.com/libprima/prima/blob/main/fortran/lincoa/trustregion.f90
    (an improvement of the original code from M. J. D. Powell's LINCOA software [Pow2015])
    
    [Pow2015] M. J. D. Powell, On Fast Trust Region Methods for Quadratic Models with Linear Constraints,
              Mathematical Programming Computation 7:3 (2015), pp. 237-267.
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H, A=A, b=b, xopt=xopt)
    if tol <= 0.0:
        raise RuntimeError("tol must be strictly positive")

    s = py_trslin(delta, g, H, A, b, xopt, tol)
    return s


def solve(delta, g, H, sl=None, su=None, A=None, b=None, xopt=None, tol=1e-2, solve_global=False, subproblem='tr'):
    """
    General interface for solving trust-region and cubic regularization subproblems:
        s = trustregion.solve(delta, g, H,
                              sl=None, su=None, xopt=None, A=None, b=None
                              tol=1e-2, solve_global=False, subproblem='tr')

    If solve_trs=True, this solves the trust-region subproblem with optional bound or linear
    constraints (but not both):
        min_{s in R^n} g^T * s + 0.5 * s^T * H * s
        s.t.  ||s||_2 <= delta
              Optional: sl <= xopt + s <= su   --OR--   A @ (xopt + s) <= b
    
    Note: the point xopt must be chosen so that s=0 is feasible
    (i.e. sl <= xopt <= su or A @ xopt <= b must hold)
    
    If solve_global=False or bound/linear constraints are provided, the function returns an approximate
    solution with the termination tolerance 'tol'. If solve_global=True and no bound/linear constraints are
    provided, a global minimizer is calculated.
    
    If subproblem=='cr', this globally solves the unconstrained cubic regularization problem:
        min_{s in R^n} g^T * s + 0.5 * s^T * H * s + (delta/3) * ||s||_2^3
    
    The global solvers are best suited to small/medium scale problems (e.g. n <= 500).

    Inputs are:
    - delta = non-negative trust-region radius
    - g = np.ndarray of shape (n,)
    - H = symmetric np.ndarray of shape (n,n)
    - sl = lower bounds, np.ndarray of shape (n,)
    - su = upper bounds, np.ndarray of shape (n,)
    - A = constraint matrix, np.ndarray of shape (m,n)
    - b = constraint RHS vector, np.ndarray of shape (m,)
    - xopt = base point for bound constraints, np.ndarray of shape (n,)
    - tol = strictly positive termination tolerance
    - solve_global = boolean, True if a global minimizer is desired
    - subproblem = 'tr' or 'cr' for trust-region or cubic regularization subproblem

    Output is:
    - s = approximate minimizer, np.ndarray of shape (n,)
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H, sl=sl, su=su, A=A, b=b, xopt=xopt)
    if tol <= 0.0:
            raise RuntimeError("tol must be strictly positive")

    if subproblem.lower() == 'tr':
        solve_trs = True
    elif subproblem.lower() == 'cr':
        solve_trs = False
    else:
        raise RuntimeError("Unknown subproblem type '%s'" % subproblem)

    if cons_type == 'none':
        if solve_trs and solve_global:
            s, lda = trsunc(delta, g, H)
        elif solve_trs:
            s, crvmin, info = trsapp(delta, g, H)
        else:
            s = arcunc(delta, g, H)

    elif cons_type == 'bounds':
        if solve_trs == False:
            raise RuntimeError("Cannot solve bound-constrained cubic regularization subproblem")
        s, crvmin = trsbox(delta, g, H, sl=sl, su=su, xopt=xopt, tol=tol)
    
    elif cons_type == 'linear':
        if solve_trs == False:
            raise RuntimeError("Cannot solve linearly constrained cubic regularization subproblem")
        s = trslin(delta, g, H, A=A, b=b, xopt=xopt, tol=tol)
    
    else:
        raise RuntimeError("Could not determine constraint type")

    return s