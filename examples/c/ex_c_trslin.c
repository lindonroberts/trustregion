#include "stdio.h"
#include "stdlib.h"
#include "trustregion.h"

int main()
{
    /*
    * TRSLIN (linearly-constrained trust-region subproblem using CG) example
    */

    const int n = 3;
    const int m = 6;
    double delta;
    double crvmin;
    double tol;
    double* g;
    double* hess;
    double* s;
    double* amat;
    double* bvec;
    double* sl;
    double* su;
    double* xbase;
    int i;

    g = (double*)malloc(sizeof(double) * n);
    hess = (double*)malloc(sizeof(double) * n * n);
    s = (double*)malloc(sizeof(double) * n);
    amat = (double*)malloc(sizeof(double) * m * n);
    bvec = (double*)malloc(sizeof(double) * m);
    sl = (double*)malloc(sizeof(double) * n);
    su = (double*)malloc(sizeof(double) * n);
    xbase = (double*)malloc(sizeof(double) * n);

    g[0] = 1.0; g[1] = 0.0; g[2] = 1.0;

    hess[0] = 1.0; hess[1] = 0.0; hess[2] = 0.0;
    hess[3] = 0.0; hess[4] = 2.0; hess[5] = 0.0;
    hess[6] = 0.0; hess[7] = 0.0; hess[8] = 2.0;

    xbase[0] = 1.0; xbase[1] = 1.0; xbase[2] = 1.0;
    sl[0] = -0.5 + xbase[0]; sl[1] = -10.0 + xbase[1]; sl[2] = -10.0 + xbase[2];
    su[0] = 11.0 + xbase[0]; su[1] = 11.0 + xbase[1]; su[2] = 11.0 + xbase[2];

    delta = 2.0;
    tol = 1.0e-2;

    // Initialize amat = 0
    for (i = 0; i < m * n; ++i) amat[i] = 0.0;

    // Change bounds sl <= xbase + s <= su 
    // into constraints amat * (xbase + s) <= bvec
    for (i = 0; i < n; ++i)
    {
        amat[i * n + i] = 1.0;
        bvec[i] = su[i];
        amat[(n + i) * n + i] = -1.0;
        bvec[n + i] = -sl[i];
    }

    f_trslin(n, delta, g, hess, m, amat, bvec, xbase, tol, s);

    printf("TRSLIN calculated original step\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);

    // Replace final constraint xbase[3] + s[3] >= -10 (unused) with a constraint forcing the 
    // bound constraint solution s=[-0.5, 0, -0.5] to not be valid
    // Use -e * (xopt + s) <= -3, so since e * xopt = 3 this means s[1] + s[2] + s[3] >= 0
    for (i = 0; i < n; ++i) amat[(m-1) * n + i] = -1.0;
    bvec[m-1] = -3.0;

    f_trslin(n, delta, g, hess, m, amat, bvec, xbase, tol, s);

    printf("TRSLIN calculated modified step\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);

    free(g);
    free(hess);
    free(s);
    free(amat);
    free(bvec);
    free(sl);
    free(su);
    free(xbase);

    return 0;
};
