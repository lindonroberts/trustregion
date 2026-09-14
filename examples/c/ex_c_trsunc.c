#include "stdio.h"
#include "stdlib.h"
#include "trustregion.h"

int main()
{
    /*
    * TRSAPP (unconstrained trust-region subproblem using CG) example
    * Example comes from [GRT2010] 
    * 
    * [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
    * regularised subproblems in optimization. Mathematical Programming Computation 2:1
    * (2010), pp. 21-57. 
    */

    const int n = 3;
    double delta;
    double lambda;
    double* g;
    double* hess;
    double* s;

    g = (double*)malloc(sizeof(double) * n);
    hess = (double*)malloc(sizeof(double) * n * n);
    s = (double*)malloc(sizeof(double) * n);

    g[0] = 5.0; g[1] = 0.0; g[2] = 4.0;

    hess[0] = 1.0; hess[1] = 0.0; hess[2] = 4.0;
    hess[3] = 0.0; hess[4] = 2.0; hess[5] = 0.0;
    hess[6] = 4.0; hess[7] = 0.0; hess[8] = 3.0;

    delta = 1.0;

    f_trsunc(n, &delta, g, hess, &lambda, s);

    printf("TRSUNC easy case\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);
    printf("lambda = %e\n", lambda);
    printf("Note: global min is s = (-1, 0, 0) with lambda=4\n");

    g[0] = 0.0; g[1] = 2.0; g[2] = 0.0;
    f_trsunc(n, &delta, g, hess, &lambda, s);

    printf("TRSUNC hard case\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);
    printf("lambda = %e\n", lambda);
    printf("Note: global min is s = (+/-0.68926566, -0.48507125, -/+0.53816237) with lambda=2.1231056\n");

    g[0] = 0.0; g[1] = 2.0; g[2] = 0.0001;
    f_trsunc(n, &delta, g, hess, &lambda, s);
    
    printf("TRSUNC nearly hard case\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);
    printf("lambda = %e\n", lambda);
    printf("Note: global min is s = (0.6892634, -0.48506297, -0.53817273) with lambda=2.1231760\n");

    free(g);
    free(hess);
    free(s);

    return 0;
};
