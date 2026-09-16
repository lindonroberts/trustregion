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
    Global unconstrained trust-region subproblem solver

    s, lambda = trustregion.trsunc(delta, g, H)
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H)
    s, lda = py_trsunc(delta, g, H)
    return s, lda


def arcunc(delta, g, H):
    """
    Global cubic regularization subproblem solver

    s = trustregion.arcunc(delta, g, H)
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H)
    s = py_arcunc(delta, g, H)
    return s

def trsapp(delta, g, H, tol=1e-2):
    """
    Approximate unconstrained trust-region subproblem solver

    s, crvmin, info = trustregion.trsapp(delta, g, H, tol=1e-2)
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H)
    if tol <= 0.0:
        raise RuntimeError("tol must be strictly positive")

    s, crvmin, info = py_trsapp(delta, g, H, tol)
    return s, crvmin, info

def trsbox(delta, g, H, sl, su, xopt, tol=1e-2):
    """
    Approximate bound-constrained trust-region subproblem solver
    
    s, crvmin = trustregion.trsbox(delta, g, H, sl, su, xopt, tol=1e-2)
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H, sl=sl, su=su, xopt=xopt)
    if tol <= 0.0:
        raise RuntimeError("tol must be strictly positive")

    s, crvmin = py_trsbox(delta, g, H, sl, su, xopt, tol)
    return s, crvmin

def trslin(delta, g, H, A, b, xopt, tol=1e-2):
    """
    Approximate linearly constrained trust-region subproblem solver
    
    s = trustregion.trslin(delta, g, H, A, b, xopt, tol=1e-2)
    """
    delta, g, H, sl, su, xopt, A, b, cons_type = validate_inputs(delta, g, H, A=A, b=b, xopt=xopt)
    if tol <= 0.0:
        raise RuntimeError("tol must be strictly positive")

    s = py_trslin(delta, g, H, A, b, xopt, tol)
    return s
