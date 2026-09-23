# trustregion: Subproblem solvers for nonlinear optimization

`trustregion` is a Fortran package (with C and Python interfaces) for solving trust region and cubic regularization subproblems.
These are common subproblems that arise in nonlinear optimization.
It is a newer version of the [trust-region](https://github.com/lindonroberts/trust-region/) Python package.

The routines `trsapp`, `trsbox` and `trslin` are taken directly from [PRIMA](https://github.com/libprima/prima), which is made available under the same BSD-3 licence as the `trustregion` package.

For a mathematical background on trust region and cubic regularization methods, see references (CGT2000), (NW2006), or (CGT2022) below.

## Installation

The main Fortran library and accompanying C interface can be built using [cmake](https://cmake.org/):

    git clone https://github.com/lindonroberts/trustregion.git
    cd trustregion
    cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/install/dir
    cmake --build build --target install

The `CMAKE_INSTALL_PREFIX` option above ensures that the headers and library will be saved in `/path/to/install/dir/include` and `/path/to/install/dir/lib` respectively.

The Python package can be installed using [pip](https://packaging.python.org/en/latest/tutorials/installing-packages/):

    git clone https://github.com/lindonroberts/trustregion.git
    cd trustregion
    pip install .

## Available solvers

**Unconstrained Trust Region Subproblem (trsunc, trsapp)**

The routines `trsunc` and `trsapp` solve the unconstrained trust region subproblem

$$\min_{s\in\mathbb{R}^n} g^T s + \frac{1}{2} s^T H s, \qquad \text{subject to} \quad \|s\|_2 \leq \Delta, $$

for vector $g\in\mathbb{R}^n$, symmetric $n\times n$ matrix $H$ and radius $\Delta \geq 0$.
The solver `trsunc` finds a global minimizer and `trsapp` finds an approximate minimizer (but is faster and better-suited to high-dimensional problems).

The global solver `trsunc` implements the method from (GRT2010) and also returns a scalar $\lambda\geq 0$ such that $(H + \lambda I)s = -g$ (i.e. the Lagrange multiplier associated with the constraint $\|s\|_2^2 \leq \Delta^2$).
The approximate solver `trsapp` implements the conjugate gradient-based method from (Pow2006).

**Unconstrained Cubic Regularization Subproblem (arcunc)**

The routine `arcunc` solves the cubic regularization subproblem

$$\min_{s\in\mathbb{R}^n} g^T s + \frac{1}{2} s^T H s +  \frac{\Delta}{3} ||s||_2^3, $$

for vector $g\in\mathbb{R}^n$, symmetric $n\times n$ matrix $H$ and radius $\Delta \geq 0$.
It finds a global minimizer using the method from (GRT2010).

**Bound Constrained Trust Region Subproblem (trsbox)**

The routine `trsbox` solves the bound-constrained trust region subproblem

$$\min_{s\in\mathbb{R}^n} g^T s + \frac{1}{2} s^T H s, \qquad \text{subject to} \quad \|s\|_2 \leq \Delta, \: \text{and} \: s_l \leq x_{opt} + s \leq s_u, $$

for vectors $g,s_l,x_{opt},s_u\in\mathbb{R}^n$, symmetric $n\times n$ matrix $H$ and radius $\Delta \geq 0$.
The solver assumes that the point $s=0$ is feasible for the bound constraints (i.e. $s_l \leq x_{opt} \leq s_u$ holds).
It implements the active set conjugate gradient method from (Pow2009).

**Linearly Constrained Trust Region Subproblem (trslin)**

The routine `trslin` solves the linearly constrained trust region subproblem (with $m$ linear inequality constraints)

$$\min_{s\in\mathbb{R}^n} g^T s + \frac{1}{2} s^T H s, \qquad \text{subject to} \quad \|s\|_2 \leq \Delta, \: \text{and} \: A (x_{opt} + s) \leq b, $$

for vectors $g,x_{opt},s_u\in\mathbb{R}^n$ and $b\in\mathbb{R}^m$, symmetric $n\times n$ matrix $H$, $m\times n$ matrix $A$ and radius $\Delta \geq 0$.
The solver assumes that the point $s=0$ is feasible for the bound constraints (i.e. $A x_{opt} \leq b$ holds).
It implements the active set conjugate gradient method from (Pow2015).

## Fortran/C usage

Fortran and C examples for each of the 5 solvers are available in the `examples/fortran` and `examples/c` directory.
These are automatically built and the executables are available in `build/examples/fortran/*` and `build/examples/c/*`, unless you set the cmake flag `BUILD_TESTING=OFF`.

The call definitions for the 5 solvers are given in `src/trustregion.F90` (Fortran) and `interfaces/c/trustregion.h` (C).

## Python usage

Python examples for each of the 5 solvers are available in the `examples/python` directory.
The full interfaces for the solvers are:

    # Unconstrained trust region subproblem
    s = trustregion.arcunc(delta, g, H)
    s, crvmin, info = trustregion.trsapp(delta, g, H, tol=1e-2)

    # Cubic regularization subproblem
    s = trustregion.arcunc(delta, g, H)

    # Bound constrained trust region subproblem
    s, crvmin = trustregion.trsbox(delta, g, H, sl, su, xopt, tol=1e-2)

    # Linearly constrained trust region subproblem
    s = trustregion.trslin(delta, g, H, A, b, xopt, tol=1e-2)

where all vector/matrix inputs are NumPy arrays (of compatible dimensions), and `tol` is a strictly positive termination tolerance for the conjugate gradient solvers.

The output `crvmin` (type `float`) from `trsapp` and `trsbox` is an estimate of the smallest eigenvalue of `H` based on the observed sequence of iterates.
The output `info` from `trsapp` is a termination flag (`info == 0` is standard termination, `info > 0` for termination due to slow progress, and `info < 0` for termination due to rounding errors).

## References

* (CGT2000) A. R. Conn, N. I. M. Gould, P. L. Toint. Trust-Region Methods. SIAM (2000).
* (CGT2022) C. Cartis, N. I. M. Gould, P. L. Toint. Evaluation Complexity of Algorithms for Nonconvex Optimization. SIAM (2022).
* (GRT2010) N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other regularised subproblems in optimization. Mathematical Programming Computation 2:1 (2010), pp. 21-57. 
* (NW2006) J. Nocedal, S. W. Wright. Nonlinear Optimization, 2nd ed. Springer (2006).
* (Pow2006) M. J. D. Powell, The NEWUOA Software for Unconstrained Optimization without Derivatives, in: Large-Scale Nonlinear Optimization, P. Pardalos, G. Di Pillo, M. Roma, eds, Springer, 2006, pp. 225-297.
* (Pow2009) M. J. D. Powell, The BOBYQA Algorithm for Bound Constrained Optimization without Derivatives, Technical report DAMTP 2009/NA06, University of Cambridge, 2006.
* (Pow2015) M. J. D. Powell, On Fast Trust Region Methods for Quadratic Models with Linear Constraints, Mathematical Programming Computation 7:3 (2015), pp. 237-267.
