!--------------------------------------------------------------------------------------------------!
! TRSAPP (unconstrained trust-region subproblem using CG) example
! Example comes from [GRT2010] 
! 
! [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other 
!           regularised subproblems in optimization. Mathematical Programming Computation 2:1
!           (2010), pp. 21-57. 
!--------------------------------------------------------------------------------------------------!


program ex_trsapp

use, intrinsic :: iso_fortran_env, only : RP => REAL64
!use consts_mod, only : RP
use trsapp_mod, only : trsapp

implicit none

integer, parameter :: n = 3
integer :: info
real(RP) :: delta, crvmin, tol, g(n), hess(n, n), s(n)

g = (/ 5.0, 0.0, 4.0 /)
hess = reshape([1.0, 0.0, 4.0, 0.0, 2.0, 0.0, 4.0, 0.0, 3.0], [n,n])
delta = 1.0
tol = 1.0E-2

call trsapp(delta, g, hess, tol, crvmin, s, info)

print *, "TRSAPP calculated step s ="
print *, s
print *, "crvmin = ", crvmin
print *, "info =", info
print *, "Note: global min is s=(-1, 0, 0)"

end program ex_trsapp