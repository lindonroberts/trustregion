!--------------------------------------------------------------------------------------------------!
! TRSLIN (linearly-constrained trust-region subproblem using CG) example
!--------------------------------------------------------------------------------------------------!


program ex_trslin

use, intrinsic :: iso_fortran_env, only : RP => REAL64
use trustregion_mod, only : trslin

implicit none

integer, parameter :: n = 3
integer, parameter :: m = 2 * n
integer :: i
real(RP) :: delta, tol
real(RP) :: g(n), hess(n, n), s(n), sl(n), su(n), xbase(n)
real(RP) :: amat(n, m), bvec(m)

g = (/ 1.0, 0.0, 1.0 /)
hess = reshape([1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0], [n,n])
xbase = 1.0
sl = (/ -0.5, -10.0, -10.0 /) + xbase
su = (/ 11.0, 11.0, 11.0 /) + xbase
delta = 2.0
tol = 1.0E-2

! Change bounds sl <= xbase + s <= su 
! into constraints amat^T * (xbase + s) <= bvec
amat = 0.0
do i = 1, n
    amat(i, i) = 1.0
    bvec(i) = su(i)
    amat(i, 2*i) = -1.0
    bvec(2*i) = -sl(i)
end do

call trslin(amat, bvec, xbase, delta, g, hess, tol, s)

print *, "TRSLIN calculated step s ="
print *, s

end program ex_trslin