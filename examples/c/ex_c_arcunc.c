#include "stdio.h"
#include "stdlib.h"
#include "trustregion.h"

int main()
{
    /*
    * ARCUNC (global unconstrained cubic regularization subproblem) example
    * Example comes from [GRT2010] 
    * 
    * [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
    * regularised subproblems in optimization. Mathematical Programming Computation 2:1
    * (2010), pp. 21-57. 
    */

    const int n = 3;
    double delta;
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

    f_arcunc(n, delta, g, hess, s);

    printf("ARCUNC easy case\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);
    printf("Note: global min is s = (-2.4827526, 0, 1.0418972)\n");  // lambda=2.6925100

    g[0] = 0.0; g[1] = 2.0; g[2] = 0.0;
    f_arcunc(n, delta, g, hess, s);

    printf("ARCUNC hard case\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);
    printf("Note: global min is s = (+/-1.62917555, -0.48507125, -/+1.27203396)\n");  //  lambda=2.1231056

    g[0] = 0.0; g[1] = 2.0; g[2] = 0.0001;
    f_arcunc(n, delta, g, hess, s);
    
    printf("ARCUNC nearly hard case\n");
    printf("s = [%e, %e, %e]\n", s[0], s[1], s[2]);
    printf("Note: global min is s = (1.62920031, -0.48506775, -1.27205329)\n"); // lambda=2.1231354

    free(g);
    free(hess);
    free(s);

    return 0;
};
