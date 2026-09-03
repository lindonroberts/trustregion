submodule (trustregion_mod) qdec_mod
!--------------------------------------------------------------------------------------------------!
! Compute the model decrease from a given step
!--------------------------------------------------------------------------------------------------!

implicit none

contains


function qdec(g, hess, s) result(dec)
!--------------------------------------------------------------------------------------------------!
! This function calculates the model decrease for a step s, given gradient G and hessian HESS
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod, inprod, issymmetric

implicit none

! Inputs
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: hess(:, :)  ! HESS(N, N)
real(RP), intent(in) :: s(:)  ! S(N)
real(RP) :: dec

! Local variables
character(len=*), parameter :: srname = 'QDEC'
integer(IK) :: n
real(RP) :: hs(size(g))

! Sizes.
n = int(size(g), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(g) == n, 'SIZE(G) == N', srname)
    call assert(size(hess, 1) == n .and. issymmetric(hess), 'HESS is n-by-n and symmetric', srname)
end if

hs = matprod(hess, s)
dec = -inprod(g, s) - HALF * inprod(s, hs)

end function qdec

function crdec(g, hess, s, delta) result(dec)
!--------------------------------------------------------------------------------------------------!
! This function calculates the cubic regularization model decrease
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP

implicit none

! Inputs
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: hess(:, :)  ! HESS(N, N)
real(RP), intent(in) :: s(:)  ! S(N)
real (RP), intent(in) :: delta
real(RP) :: dec

dec = qdec(g, hess, s) - (delta / 3.0) * sqrt(sum(s**2))**3

end function qdec


end submodule qdec_mod