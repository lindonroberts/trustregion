#include "stdio.h"
#include "stdlib.h"
#include "trustregion.h"

int main()
{
    /*
    * TRSBOX (box-constrained trust-region subproblem using CG) example
    * Example comes from [GRT2010] 
    * 
    * [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
    * regularised subproblems in optimization. Mathematical Programming Computation 2:1
    * (2010), pp. 21-57. 
    */

    const int n = 3;
    double delta;
    double crvmin;
    double tol;
    double* g;
    double* hess;
    double* s;
    double* sl;
    double* su;
    double* xbase;

    g = (double*)malloc(sizeof(double) * n);
    hess = (double*)malloc(sizeof(double) * n * n);
    s = (double*)malloc(sizeof(double) * n);
    sl = (double*)malloc(sizeof(double) * n);
    su = (double*)malloc(sizeof(double) * n);
    xbase = (double*)malloc(sizeof(double) * n);

    g[0] = 5.0;
    g[1] = 0.0;
    g[2] = 4.0;

    hess[0] = 1.0; hess[1] = 0.0; hess[2] = 4.0;
    hess[3] = 0.0; hess[4] = 2.0; hess[5] = 0.0;
    hess[6] = 4.0; hess[7] = 0.0; hess[8] = 3.0;

    xbase[0] = 0.0; xbase[1] = 0.0; xbase[2] = 0.0;
    sl[0] = -0.1; sl[1] = -1.0; sl[2] = -1.0;
    su[0] = 1.0; su[1] = 1.0; su[2] = 1.0;

    delta = 1.0;
    tol = 1.0e-2;

    f_trsbox(n, &delta, g, hess, sl, su, &tol, xbase, &crvmin, s);

    printf("TRSBOX calculated step\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);
    printf("crvmin = %e\n", crvmin);
    printf("Note: global min (no box) is s = (-1, 0, 0)\n");

    free(g);
    free(hess);
    free(s);
    free(sl);
    free(su);
    free(xbase);

    return 0;
};
