!--------------------------------------------------------------------------------------------------!
! TRSBOX (box-constrained trust-region subproblem using CG) example
! Example comes from [GRT2010] 
! 
! [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
!           regularised subproblems in optimization. Mathematical Programming Computation 2:1
!           (2010), pp. 21-57. 
!--------------------------------------------------------------------------------------------------!


program ex_trsapp

use, intrinsic :: iso_fortran_env, only : RP => REAL64
use trustregion_mod, only : trsbox

implicit none

integer, parameter :: n = 3
real(RP) :: delta, crvmin, tol, g(n), hess(n, n), s(n), sl(n), su(n), xbase(n)

g = (/ 5.0, 0.0, 4.0 /)
hess = reshape([1.0, 0.0, 4.0, 0.0, 2.0, 0.0, 4.0, 0.0, 3.0], [n,n])
xbase = 0.0
sl = (/ -0.1, -1.0, -1.0 /)
su = (/ 1.0, 1.0, 1.0 /)
delta = 1.0
tol = 1.0E-2

call trsbox(delta, g, hess, sl, su, tol, xbase, crvmin, s)

print *, "TRSBOX calculated step s ="
print *, s
print *, "crvmin = ", crvmin
print *, "Note: global min (no box) is s=(-1, 0, 0)"

end program ex_trsapp