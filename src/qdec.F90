module qdec_mod
!--------------------------------------------------------------------------------------------------!
! Compute the model decrease from a given step
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: quadform, qdec, crdec

contains

function quadform(hess, s) result(q)
    !------------------------!
    ! Compute dot(s, hess*s)
    !------------------------!

    ! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : matprod, inprod, issymmetric

    implicit none

    ! Inputs
    real(RP), intent(in) :: hess(:, :)  ! HESS(N, N)
    real(RP), intent(in) :: s(:)  ! S(N)
    real(RP) :: q

    ! Local variables
    character(len=*), parameter :: srname = 'QDEC'
    integer(IK) :: n
    integer(IK) :: i, j

    ! Sizes.
    n = int(size(hess, 1), kind(n))

    ! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(s) == n, 'SIZE(S) == N', srname)
        call assert(size(hess, 1) == n .and. issymmetric(hess), 'HESS is n-by-n and symmetric', srname)
    end if

    q = ZERO
    do i = 1, n 
        do j = 1, n 
            q = q + hess(i, j) * s(i) * s(j)
        end do
    end do
end function quadform

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

end function crdec


end module qdec_mod