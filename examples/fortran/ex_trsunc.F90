!--------------------------------------------------------------------------------------------------!
! TRSUNC (global unconstrained trust-region subproblem) example
! Example comes from [GRT2010] 
! 
! [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
!           regularised subproblems in optimization. Mathematical Programming Computation 2:1
!           (2010), pp. 21-57. 
!--------------------------------------------------------------------------------------------------!


program ex_trsapp

use, intrinsic :: iso_fortran_env, only : RP => REAL64
use trustregion_mod, only : trsunc

implicit none

integer, parameter :: n = 3
integer :: info
real(RP) :: delta, lambda, g(n), hess(n, n), s(n), strue(n), strue2(n), lambda_true

! Easy case
g = (/ 5.0, 0.0, 4.0 /)
hess = reshape([1.0, 0.0, 4.0, 0.0, 2.0, 0.0, 4.0, 0.0, 3.0], [n,n])
delta = 1.0

strue = (/ -1, 0, 0 /)
lambda_true = 4

call trsunc(delta, g, hess, lambda, s)

print *, "TRSUNC easy case"
print *, "s = ", s
print *, "lambda = ", lambda
print *, "strue =", strue
print *, "lambda_true = ", lambda_true

! Hard case
g = (/ 0.0, 2.0, 0.0 /)
delta = 1.0
strue = (/ -0.68926566, -0.48507125, 0.53816237 /)
strue2 = (/ 0.68926566, -0.48507125, -0.53816237 /)
lambda_true = 2.1231056256176606

call trsunc(delta, g, hess, lambda, s)

print *, ""
print *, "TRSUNC hard case"
print *, "s = ", s
print *, "lambda = ", lambda
print *, "strue =", strue
print *, "strue2 =", strue2
print *, "lambda_true = ", lambda_true

! Almost hard case

g = (/ 0.0, 2.0, 0.0001 /)
delta = 1.0
strue = (/ 0.6892634, -0.48506297, -0.53817273 /)
lambda_true = 2.123176000326642

call trsunc(delta, g, hess, lambda, s)

print *, ""
print *, "TRSUNC almost hard case"
print *, "s = ", s
print *, "lambda = ", lambda
print *, "strue =", strue
print *, "lambda_true = ", lambda_true

end program ex_trsapp