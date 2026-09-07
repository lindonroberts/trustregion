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
real(RP) :: delta, lambda, tol, g(n), hess(n, n), s(n)

g = (/ 5.0, 0.0, 4.0 /)
hess = reshape([1.0, 0.0, 4.0, 0.0, 2.0, 0.0, 4.0, 0.0, 3.0], [n,n])
delta = 1.0
tol = 1.0E-2

call trsunc(delta, g, hess, lambda, s)

print *, "TRSUNC calculated step s ="
print *, s
print *, "lambda = ", lambda
print *, "Note: global min is s=(-1, 0, 0)"

end program ex_trsapp