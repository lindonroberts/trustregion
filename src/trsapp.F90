module trsapp_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of NEWUOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Saturday, April 06, 2024 PM11:03:07
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trsapp

contains


subroutine trsapp(delta, g_in, hess_in, tol, crvmin, s, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds an approximate solution to the N-dimensional trust region subproblem
!
! min <S, GOPT> + 0.5*<S, HESS*S> s.t. ||S|| <= DELTA
!
! The calculation of S begins with the truncated conjugate gradient method. If the boundary of the
! trust region is reached, then further changes to S may be made, each one being in the 2-dimensional
! space spanned by the current S and the corresponding gradient of Q. Thus S should provide a
! substantial reduction to Q within the trust region. See Section 5 of the NEWUOA paper.
!
! At return, S will be the approximate solution. CRVMIN will be set to the least curvature of
! HESS along the conjugate directions that occur, except that it is set to ZERO if S goes all the
! way to the trust-region boundary. INFO is an exit flag with the following possible values.
! - INFO = 0: an approximate solution satisfying one of the following conditions is found:
! 1. ||G+HS||/||G0|| <= TOL,
! 2. ||S|| = DELTA and <S, -(G+HS)> >= (1 - TOL)*||S||*||G+HS||,
! where TOL is set to 1e-2 in NEWUOA;
! - INFO = 1: the last iteration reduces Q only insignificantly;
! - INFO = 2: the maximal number of iterations is attained;
! - INFO = -1: too much rounding error to continue.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, HALF, REALMIN, ZERO, TENTH, EPS, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : inprod, issymmetric, norm, project, matprod
use, non_intrinsic :: univar_mod, only : circle_min

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: g_in(:)   ! G_IN(N)
real(RP), intent(in) :: hess_in(:, :)    ! HESS_IN(N, N)
real(RP), intent(in) :: tol

! Outputs
real(RP), intent(out) :: crvmin
real(RP), intent(out) :: s(:)   ! S(N)
integer(IK), intent(out), optional :: info

! Local variables
character(len=*), parameter :: srname = 'TRSAPP'
integer(IK) :: info_loc
integer(IK) :: iter
integer(IK) :: maxiter
integer(IK) :: n
logical :: scaled
logical :: twod_search
real(RP) :: alpha
real(RP) :: angle
real(RP) :: args(4)
real(RP) :: bstep
real(RP) :: cth
real(RP) :: d(size(g_in))
real(RP) :: dd
real(RP) :: delsq
real(RP) :: dg
real(RP) :: dhd
real(RP) :: dhs
real(RP) :: ds
real(RP) :: g(size(g_in))
real(RP) :: gg
real(RP) :: gg0
real(RP) :: ggsav
real(RP) :: gopt(size(g_in))
real(RP) :: hd(size(g_in))
real(RP) :: hess(size(hess_in, 1), size(hess_in, 2))
real(RP) :: hs(size(g_in))
real(RP) :: modscal
real(RP) :: qadd
real(RP) :: qred
real(RP) :: reduc
real(RP) :: resid
real(RP) :: sg
real(RP) :: shs
real(RP) :: sold(size(g_in))
real(RP) :: sqrtd
real(RP) :: ss
real(RP) :: sth

! Sizes
n = int(size(g_in), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(size(g_in) == n, 'SIZE(G) = N', srname)
    call assert(size(hess_in, 1) == n .and. issymmetric(hess_in), 'HESS is an NxN symmetric matrix', srname)
    call assert(size(s) == n, 'SIZE(S) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

! Scale the problem if GOPT contains large values. Otherwise, floating point exceptions may occur.
! Note that CRVMIN must be scaled back if it is nonzero, but the step is scale invariant.
! N.B.: It is faster and safer to scale by multiplying a reciprocal than by division. See
! https://fortran-lang.discourse.group/t/ifort-ifort-2021-8-0-1-0e-37-1-0e-38-0/
if (maxval(abs(g_in)) > 1.0E12) then   ! The threshold is empirical.
    modscal = max(TWO * REALMIN, ONE / maxval(abs(g_in)))  ! MAX: precaution against underflow.
    gopt = g_in * modscal
    hess = hess_in * modscal
    scaled = .true.
else
    modscal = ONE  ! This value is not used, but Fortran compilers may complain without it.
    gopt = g_in
    hess = hess_in
    scaled = .false.
end if

s = ZERO
crvmin = ZERO
qred = ZERO
info_loc = 2  ! Default exit flag is 2, i.e., MAXITER is attained

! Prepare for the first line search.
!--------------------------------------------------------------------------------------------------!
! N.B.: During the iterations, G is NOT updated, and it equals always GOPT, which is the gradient
! of the trust-region model at the trust-region center X. However, GG is updated: GG = ||G + HS||^2,
! which is the norm square of the gradient at the current iterate.
g = gopt
gg = inprod(g, g)
!--------------------------------------------------------------------------------------------------!
gg0 = gg
d = -g
dd = gg
ds = ZERO
ss = ZERO
hs = ZERO
delsq = delta * delta
maxiter = n

twod_search = .false.

! The truncated-CG iterations.
!
! The iteration will be terminated in 4 possible cases:
! 1. the maximal number of iterations is attained;
! 2. QADD <= TOL*QRED or ||G|| <= TOL*||G0||, where QADD is the reduction of Q due to the latest
! CG step, QRED is the reduction of Q since the beginning until the latest CG step, G is the
! current gradient, and G0 is the initial gradient; see (5.13) of the NEWUOA paper;
! 3. DS <= 0
! 4. ||S|| = DELTA, i.e., CG path cuts the trust region boundary.
!
! In the 4th case, twod_search will be set to true, meaning that S will be improved by a sequence of
! two-dimensional search, the two-dimensional subspace at each iteration being span(S, -G).
do iter = 1, maxiter
    ! Exit if G contains NaN.
    if (is_nan(gg)) then
        info_loc = -1
        exit
    end if
    ! Exit if GG is small. This must be done first; otherwise, DD can be 0 and BSTEP is not well
    ! defined. The inequality below must be non-strict so that GG = GG0 = 0 will trigger the exit.
    if (gg <= (tol**2) * gg0) then
        info_loc = 0
        exit
    end if

    ! Set BSTEP to the step length such that ||S + BSTEP*D|| = DELTA.
    if (iter == 1) then
        bstep = delta / sqrt(dd)
    else
        resid = delsq - ss
        if (resid <= 0) then
            twod_search = .true.
            exit
        end if

        ! Powell's code does not have the following two IFs.
        !--------------------------------------------------!
        if (dd <= EPS * delsq) then
            info_loc = 0
            exit
        end if
        if (is_nan(ds)) then
            info_loc = -1
            exit
        end if
        !--------------------------------------------------!

        ! SQRTD: square root of a discriminant. The MAXVAL avoids SQRTD < ABS(DS) due to underflow.
        sqrtd = maxval([sqrt(ds**2 + dd * resid), abs(ds), sqrt(dd * resid)])
        ! Powell's code does not distinguish the following two cases, which have no difference in
        ! precise arithmetic. The following scheme stabilizes the calculation. Copied from LINCOA.
        if (ds <= 0) then
            bstep = (sqrtd - ds) / dd
        else
            bstep = resid / (sqrtd + ds)
        end if
    end if

    ! BSTEP < 0 should not happen. BSTEP may be 0 or NaN if, e.g., DS or DD becomes Inf.
    if (bstep <= 0) then
        exit
    end if
    if (.not. is_finite(bstep)) then
        info_loc = -1
        exit
    end if

    hd = matprod(hess, d)
    dhd = inprod(d, hd)

    ! Set the step-length ALPHA and update CRVMIN.
    if (dhd <= 0) then
        alpha = bstep
    else
        alpha = min(bstep, gg / dhd)
        if (iter == 1) then
            crvmin = dhd / dd
        else
            crvmin = min(crvmin, dhd / dd)
        end if
    end if
    ! QADD is the reduction of Q due to the new CG step.
    qadd = alpha * (gg - HALF * alpha * dhd)
    ! QRED is the reduction of Q up to now.
    qred = qred + qadd
    ! QADD and QRED will be used in the 2-dimensional minimization if any.

    ! Update S, HS, and GG.
    sold = s
    s = s + alpha * d
    ss = inprod(s, s)
    hs = hs + alpha * hd
    ggsav = gg  ! Gradient norm square before this iteration
    gg = inprod(g + hs, g + hs)  ! Current gradient norm square
    ! We may record g+hs for later usage:
    ! gnew = g + hs
    ! Note that we should NOT set g = g + hs, because g contains the gradient of Q at X.

    ! Check whether to exit. This should be done after updating HS and GG, which will be used for
    ! the 2-dimensional minimization if any.
    ! Exit in case of Inf/NaN in S. This should come the first! Otherwise, we may return an S that
    ! contains NaN and fulfills other exit conditions.
    if (.not. is_finite(sum(abs(s)))) then
        s = sold
        info_loc = -1
        exit
    end if

    ! Exit if CG path cuts the boundary. It is the only possibility that TWOD_SEARCH is true.
    if (alpha >= bstep .or. ss >= delsq) then
        crvmin = ZERO
        twod_search = (n >= 2 .and. gg > (tol**2) * gg0)  ! TWOD_SEARCH should be FALSE if N = 1.
        exit
    end if

    ! Exit due to small QADD.
    if (qadd <= tol * qred) then
        info_loc = 1
        exit
    end if

    ! Prepare for the next CG iteration.
    d = (gg / ggsav) * d - g - hs  ! CG direction
    dd = inprod(d, d)
    ds = inprod(d, s)
    if (ds <= 0) then
        ! DS is positive in theory.
        info_loc = -1
        exit
    end if
end do

if (ss <= 0 .or. is_nan(ss)) then
    ! This may occur for ill-conditioned problems due to rounding.
    info_loc = -1
    twod_search = .false.
end if

if (twod_search) then
    ! At least 1 iteration of 2-dimensional minimization
    maxiter = max(1_IK, maxiter - iter)
else
    maxiter = 0
end if

! The 2-dimensional minimization
! N.B.: During the iterations, G is NOT updated, and it equals always GOPT, which is the gradient
! of the trust-region model at the trust-region center X. However, GG is updated: GG = ||G + HS||^2,
! which is the norm square of the gradient at the current iterate.
do iter = 1, maxiter
    ! Exit if G contains NaN.
    if (is_nan(gg)) then
        info_loc = -1
        exit
    end if
    ! Exit if GG is small. The inequality must be non-strict so that GG = GG0 = 0 triggers the exit.
    if (gg <= (tol**2) * gg0) then
        info_loc = 0
        exit
    end if
    sg = inprod(s, g)
    shs = inprod(s, hs)

    ! Begin the 2-dimensional minimization by calculating D and HD and some scalar products.

    ! Powell's code calculates D as follows. In precise arithmetic, INPROD(D, S) = 0, ||D|| = ||S||.
    ! However, when DELSQ*GG - SGK**2 is tiny, the error in D can be large and hence damage these
    ! equalities significantly. This did happen in tests, especially when using the single precision.
    ! !sgk = sg + shs
    ! !if (sgk / sqrt(gg * delsq) <= tol - ONE) then
    ! !    info_loc = 0
    ! !    exit
    ! !end if
    ! !t = sqrt(delsq * gg - sgk**2)
    ! !d = (delsq / t) * (g + hs) - (sgk / t) * s

    ! We calculate D as below. It did improve the performance of NEWUOA in our test.
    ! PROJECT(X, V) returns the projection of X to SPAN(V): X'*(V/||V||)*(V/||V||).
    d = (g + hs) - project(g + hs, s)
    ! N.B.:
    ! 1. The condition ||D||<=SQRT(TOL*GG) below is equivalent to |INPROD(G+HS,S)|<=SQRT((1-TOL)*GG*SS).
    ! As given above, Powell's code triggers an exit if INPROD(G+HS,S)=SGK<=(TOL-1)*SQRT(GG*SS).
    ! Since |SQRT(1-TOL) - (1-TOL)| <= TOL/2, our condition is close to |SGK| <= (TOL-1)*GG*SS.
    ! When |SGK| is tiny, S and G+HS are nearly parallel and hence the 2-dimensional search cannot
    ! continue. Note that SGK is unlikely positive if everything goes well.
    ! 2. SQRT(TOL)*SQRT(GG) is less likely to encounter underflow than SQRT(TOL*GG).
    ! 3. The condition below should be non-strict so that ||D|| = 0 can trigger the exit.
    if (norm(d) <= sqrt(tol) * sqrt(gg)) then
        info_loc = 0
        exit
    end if
    d = (norm(s) / norm(d)) * d
    ! In precise arithmetic, INPROD(D, S) = 0 and ||D|| = ||S|| = DELTA.
    if (abs(inprod(d, s)) >= TENTH * norm(d) * norm(s) .or. norm(d) >= TWO * delta) then
        info_loc = -1
        exit
    end if

    hd = matprod(hess, d)

    ! Seek the value of the angle that minimizes Q.
    ! First, calculate the coefficients of the objective function on the circle.
    dg = inprod(d, g)
    dhd = inprod(hd, d)
    dhs = inprod(hd, s)
    args = [sg, HALF * (shs - dhd), dg, dhs]
    ! The 50 in the line below was chosen by Powell. It works the best in tests, MAGICALLY. Larger
    ! (e.g., 60, 100) or smaller (e.g., 20, 40) values will worsen the performance of NEWUOA. Why??
    angle = circle_min(circle_fun_trsapp, args, 50_IK)

    ! Calculate the new S.
    cth = cos(angle)
    sth = sin(angle)
    sold = s
    s = cth * s + sth * d

    ! Exit in case of Inf/NaN in S.
    if (.not. is_finite(sum(abs(s)))) then
        s = sold
        info_loc = -1
        exit
    end if

    ! Test for convergence.
    reduc = circle_fun_trsapp(ZERO, args) - circle_fun_trsapp(angle, args)
    qred = qred + reduc
    if (reduc / qred <= tol) then
        info_loc = 1
        exit
    end if

    ! Calculate HS.
    hs = cth * hs + sth * hd
    gg = inprod(g + hs, g + hs)
end do

! Set CRVMIN to zero if it is NaN, which may happen if the problem is ill-conditioned.
if (is_nan(crvmin)) then
    crvmin = ZERO
end if

! Scale CRVMIN back before return. Note that the trust-region step is scale invariant.
if (scaled .and. crvmin > 0) then
    crvmin = crvmin / modscal
end if

if (present(info)) then
    info = info_loc
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(s) == n .and. all(is_finite(s)), 'SIZE(S) == N, S is finite', srname)
    ! Due to rounding, it may happen that ||S|| > DELTA, but ||S|| > 2*DELTA is highly improbable.
    call assert(norm(s) <= TWO * delta, '||S|| <= 2*DELTA', srname)
    call assert(crvmin >= 0, 'CRVMIN >= 0', srname)
end if

end subroutine trsapp


function circle_fun_trsapp(theta, args) result(f)
!--------------------------------------------------------------------------------------------------!
! This function defines the objective function of the 2-dimensional search on a circle in TRSAPP.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none
! Inputs
real(RP), intent(in) :: theta
real(RP), intent(in) :: args(:)

! Outputs
real(RP) :: f

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_FUN_TRSAPP'
real(RP) :: cth
real(RP) :: sth

! Preconditions
if (DEBUGGING) then
    call assert(size(args) == 4, 'SIZE(ARGS) == 4', srname)
end if

!====================!
! Calculation starts !
!====================!

cth = cos(theta)
sth = sin(theta)
f = (args(1) + args(2) * cth) * cth + (args(3) + args(4) * cth) * sth

!====================!
!  Calculation ends  !
!====================!
end function circle_fun_trsapp


end module trsapp_mod