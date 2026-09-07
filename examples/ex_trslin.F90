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
    amat(i, n+i) = -1.0
    bvec(n+i) = -sl(i)
end do

call trslin(amat, bvec, xbase, delta, g, hess, tol, s)

print *, "TRSLIN calculated original step s ="
print *, s

! Replace final constraint xbase[3] + s[3] >= -10 (unused) with a constraint forcing the 
! bound constraint solution s=[-0.5, 0, -0.5] to not be valid
! Use -e * (xopt + s) <= -3, so since e * xopt = 3 this means s[1] + s[2] + s[3] >= 0
amat(:, m) = -1.0
bvec(m) = -3.0

call trslin(amat, bvec, xbase, delta, g, hess, tol, s)

print *, "TRSLIN calculated modified step s ="
print *, s

end program ex_trslin