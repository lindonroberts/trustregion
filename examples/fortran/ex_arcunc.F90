!--------------------------------------------------------------------------------------------------!
! ARCUNC (global unconstrained cubic regularization subproblem) example
! Example comes from [GRT2010] 
! 
! [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
!           regularised subproblems in optimization. Mathematical Programming Computation 2:1
!           (2010), pp. 21-57. 
!--------------------------------------------------------------------------------------------------!


program ex_trsapp

use, intrinsic :: iso_fortran_env, only : RP => REAL64
use trustregion_mod, only : arcunc

implicit none

integer, parameter :: n = 3
integer :: info
real(RP) :: delta, lambda, g(n), hess(n, n), s(n), strue(n), strue2(n), lambda_true

! Easy case
g = (/ 5.0, 0.0, 4.0 /)
hess = reshape([1.0, 0.0, 4.0, 0.0, 2.0, 0.0, 4.0, 0.0, 3.0], [n,n])
delta = 1.0

strue = (/ -2.4827526, 0.0,  1.04189722 /)
lambda_true = 2.6925100362713921

call arcunc(delta, g, hess, lambda, s)

print *, "ARCUNC easy case"
print *, "s = ", s
print *, "lambda = ", lambda
print *, "strue =", strue
print *, "lambda_true = ", lambda_true

! Hard case
g = (/ 0.0, 2.0, 0.0 /)
delta = 1.0
strue = (/ 1.62917555, -0.48507125, -1.27203396 /)
strue2 = (/ -1.62917555, -0.48507125, 1.27203396 /)
lambda_true = 2.1231056256176606 
! alpha = 2.0669502605991603, vmin = [+/-0.78820544, 0.0, -/+0.61541221]

call arcunc(delta, g, hess, lambda, s)

print *, ""
print *, "ARCUNC hard case"
print *, "s = ", s
print *, "lambda = ", lambda
print *, "strue =", strue
print *, "strue2 =", strue2
print *, "lambda_true = ", lambda_true

! Almost hard case

g = (/ 0.0, 2.0, 0.0001 /)
delta = 1.0
strue = (/ 1.62920031, -0.48506775, -1.27205329 /)
lambda_true = 2.1231353990897641

call arcunc(delta, g, hess, lambda, s)

print *, ""
print *, "ARCUNC almost hard case"
print *, "s = ", s
print *, "lambda = ", lambda
print *, "strue =", strue
print *, "lambda_true = ", lambda_true

end program ex_trsapp