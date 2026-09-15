submodule (trustregion_mod) trslin_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of LINCOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Saturday, March 09, 2024 PM12:09:39
!--------------------------------------------------------------------------------------------------!

implicit none

contains


subroutine trslin(delta, g_in, hess_in, amat_in, bvec_in, xopt, tol, s)
!--------------------------------------------------------------------------------------------------!
! This subroutine solves
!       minimize Q(XOPT + D)  s.t. ||D|| <= DELTA, AMAT^T*(XOPT+D) <= B.
! It is assumed that D = 0 is feasible, namely AMAT^T * XOPT <= B except for rounding errors. See Powell 2015
! for details.
!
! AMAT, B, GOPT, and HESS are the same as the terms with these
! names in LINCOB.
!
! S is the total calculated step so far from the trust region centre, its final value being given by
!   the sequence of CG iterations, which terminate if the trust region boundary is reached.
! G is always the gradient of the model at the current S.
! D is the search direction of each line search.
! RESCON: If RESCON(J) is negative, then |RESCON(J)| must be no less than the trust region radius,
!   so that the J-th constraint can be ignored.
! RESNEW: A negative value of RESNEW(J) indicates that the J-th constraint does not restrict the CG
!   steps of the current trust region calculation, a zero value of RESNEW(J) indicates that the J-th
!   constraint is active, and otherwise RESNEW(J) is set to the greater of TINYCV and the actual
!   residual of the J-th constraint for the current S.
! RESACT holds the residuals of the active constraints, which may be positive.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, TWO, HALF, TEN, MAXPOW10, EPS, REALMIN, TINYCV, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm, solve, isorth, istriu, &
    & issymmetric, trueloc, matprod, eye, linspace

! Solver-specific modules
use, non_intrinsic :: getact_mod, only : getact

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: g_in(:)  ! G_IN(N)
real(RP), intent(in) :: hess_in(:, :)  ! HESS_IN(N, N)
real(RP), intent(in) :: amat_in(:, :)  ! AMAT_IN(N, M)
real(RP), intent(in) :: bvec_in(:)  ! BVEC_IN(M)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: tol

! Outputs
real(RP), intent(out) :: s(:)  ! S(N)

! Local variables
character(len=*), parameter :: srname = 'TRSLIN'
integer(IK) :: iact(size(amat_in, 2))
integer(IK) :: iter
integer(IK) :: itercg
integer(IK) :: jsav
integer(IK) :: m
integer(IK) :: maxiter
integer(IK) :: n
integer(IK) :: nact
integer(IK) :: npt
logical :: newact
real(RP) :: ad(size(amat_in, 2))
real(RP) :: alpha
real(RP) :: alphm
real(RP) :: alpht
real(RP) :: amat(size(amat_in, 1), size(amat_in, 2))
real(RP) :: anorm(size(amat_in, 2))
real(RP) :: beta
real(RP) :: bvec(size(amat_in, 2))
real(RP) :: d(size(g_in))
real(RP) :: dd
real(RP) :: delsq
real(RP) :: dg
real(RP) :: dhd
real(RP) :: dproj(size(g_in))
real(RP) :: ds
real(RP) :: frac(size(amat_in, 2))
real(RP) :: g(size(g_in))
real(RP) :: gamma
real(RP) :: gopt(size(g_in))
real(RP) :: hd(size(g_in))
real(RP) :: hess(size(hess_in, 1), size(hess_in, 2))
real(RP) :: modscal
real(RP) :: orthtol
real(RP) :: pg(size(g_in))
real(RP) :: psd(size(g_in))
real(RP) :: qfac(size(g_in), size(g_in))
real(RP) :: reduct
real(RP) :: resact(size(amat_in, 2))
real(RP) :: rescon(size(amat_in, 2))
real(RP) :: resid
real(RP) :: resnew(size(amat_in, 2))
real(RP) :: restmp(size(amat_in, 2))
real(RP) :: rfac(size(g_in), size(g_in))
real(RP) :: sold(size(s))
real(RP) :: sqrtd
real(RP) :: ss

! Sizes.
m = int(size(amat_in, 2), kind(m))
n = int(size(g_in), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(size(amat_in, 1) == n .and. size(amat_in, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(size(bvec_in) == m, 'SIZE(BVEC) == M', srname)
    call assert(size(xopt) == n, 'SIZE(XOPT) == N', srname)
    call assert(size(g_in) == n, 'SIZE(G) == N', srname)
    call assert(size(hess_in, 1) == n .and. issymmetric(hess_in), 'HESS is n-by-n and symmetric', srname)
end if

!===============!
! Preprocessing !
!===============!
! Normalize the linear constraints so that each constraint has a gradient of norm 1.
! getact assumes this
anorm = sqrt(sum(amat_in**2, dim=1))
amat = amat_in / spread(anorm, dim=1, ncopies=n)
bvec = bvec_in / anorm

! RESCON holds information about the constraint residuals at the current trust region center XOPT.
!   1. If if B(J) - AMAT(:, J)^T*XOPT <= DELTA, then RESCON(J) = B(J) - AMAT(:, J)^T*XOPT. Note that
!   RESCON >= 0 in this case, because the algorithm keeps XOPT to be feasible.
!   2. Otherwise, RESCON(J) is a negative value that B(J) - AMAT(:,J)^T*XOPT >= |RESCON(J)| >= DELTA.
!   RESCON can be updated without calculating the constraints that are far from being active, so
!   that we only need to evaluate the constraints that are nearly active.
rescon = max(bvec - matprod(xopt, amat), ZERO)
rescon(trueloc(rescon >= delta)) = -rescon(trueloc(rescon >= delta))

! Set some more initial values.
qfac = eye(n)
rfac = ZERO
nact = 0
iact = linspace(1_IK, m, m)

!====================!
! Calculation starts !
!====================!

! Scale the problem if GOPT contains large values. Otherwise, floating point exceptions may occur.
! Note that the trust-region step is scale invariant.
! N.B.: It is faster and safer to scale by multiplying a reciprocal than by division. See
! https://fortran-lang.discourse.group/t/ifort-ifort-2021-8-0-1-0e-37-1-0e-38-0/
if (maxval(abs(g_in)) > 1.0E12) then   ! The threshold is empirical.
    modscal = max(TWO * REALMIN, ONE / maxval(abs(g_in)))  ! MAX: precaution against underflow.
    gopt = g_in * modscal
    hess = hess_in * modscal
else
    gopt = g_in
    hess = hess_in
end if

! Return if G is not finite. Otherwise, GETACT will fail in the debugging mode.
if (.not. is_finite(sum(abs(gopt)))) then
    s = ZERO
    return
end if

! Set the initial elements of RESNEW, RESACT and S.

! 1. RESNEW(J) < 0 indicates that the J-th constraint does not restrict the CG steps of the current
! trust region calculation. In other words, RESCON >= DELTA.
! 2. RESNEW(J) = 0 indicates that J is an entry of IACT(1:NACT).
! 3. RESNEW(J) > 0 means that RESNEW(J) = max(B(J) - AMAT(:, J)^T*(XOPT+S), TINYCV), where S is the
! step up to now, calculated by a sequence of (truncated) CG iterations.
! N.B.: The order of the following lines is important, as the later ones override the earlier.
resnew = rescon
resnew(trueloc(rescon >= 0)) = max(TINYCV, rescon(trueloc(rescon >= 0)))
resnew(trueloc(rescon >= delta)) = -ONE
!!MATLAB:
!!resnew = rescon; resnew(rescon >= 0) = max(TINYCV, rescon(rescon >= 0)); resnew(rescon >= delta) = -1;
resnew(iact(1:nact)) = ZERO

! RESACT contains the constraint residuals of the constraints in IACT(1:NACT), namely the values
! of B(J) - AMAT(:, J)^T*(XOPT+S) for the J in IACT(1:NACT). Here, IACT(1:NACT) is a set of
! indicates such that the columns of AMAT(:, IACT(1:NACT)) form a basis of the constraint gradients
! in the "active set". For the definition of the "active set", see (3.5) of Powell (2015) and the
! comments at the beginning of the GETACT subroutine.
! N.B.: Between two calls of GETACT, S is updated in the orthogonal complement of the "active"
! gradients (i.e., null space of the "active" constraints). Therefore, RESACT remains unchanged.
! RESACT is changed right after GETACT is called if the first search direction D is not PSD
! but PSD + GAMMA * DPROJ.
resact(1:nact) = rescon(iact(1:nact))

g = gopt
delsq = delta * delta
s = ZERO
ss = ZERO
reduct = ZERO
newact = .true.

! ITERCG is the number of CG iterations corresponding to the current "active set" obtained by
! calling GETACT. These CG iterations are restricted in the orthogonal complement of the active
! gradients (i.e., null space of the active constraints).
! The following initial value of ITERCG is an artificial value that is not used. It is to entertain
! Fortran compilers (can it be be removed?).
itercg = -1

! What is the THEORETICAL upper bound of ITER? For the moment, we set the following MAXITER.
! The formulation of MAXITER below contains a precaution against overflow. In MATLAB/Python/Julia/R,
! we can write maxiter = min(10000, 10*(m + n))
maxiter = int(min(10**min(4, range(0_IK)), 10 * int(m + n)), IK)
do iter = 1, maxiter  ! Powell's code is essentially a DO WHILE loop. We impose an explicit MAXITER.
    if (newact) then
        ! GETACT picks the active set for the current S. It also sets PSD to the vector closest to
        ! -G that is orthogonal to the normals of the active constraints. PSD is scaled to have
        ! length 0.2*DELTA. Then a move of PSD from S is allowed by the linear constraints: PSD
        ! reduces the values of the nearly active constraints; it changes the inactive constraints
        ! by at most 0.2*DELTA, but the residuals of these constraints at no less than 0.2*DELTA.
        ! N.B.: The magic number 0.2 appears also in GETACT (TDEL = 0.2_RP * DELTA). It works well.
        call getact(amat, delta, g, iact, nact, qfac, resact, resnew, rfac, psd)
        dd = inprod(psd, psd)
        if (dd <= EPS * delsq .or. is_nan(dd)) then  ! Powell's code: IF (DD <= 0) THEN
            exit
        end if
        psd = (0.2_RP * delta / sqrt(dd)) * psd

        ! If the modulus of the residual of an "active constraint" is substantial (i.e., more than
        ! 1.0E-4*DELTA), then modify the searching direction PSD by a projection step to the
        ! boundaries of the "active constraint". This modified step will reduce the constraint
        ! residuals of the "active constraints" (see the update of RESACT below). The motivation is
        ! that the constraints in the "active set" are presumed to be active, and hence should have
        ! zero residuals (no constraint is violated, as the current method is feasible). According
        ! to a test on 20220821, this modification is important for the performance of LINCOA.
        ! N.B.:
        ! 1. The residual of the constraint A*X <= B is defined as B - A*X. It is not the constraint
        ! violation. Indeed, the constraint violations of the iterates are 0 in the current method.
        ! 2. We prefer `ANY(X > Y)` to `MAXVAL(X) > Y`, as Fortran standards do not specify
        ! MAXVAL(X) when X contains NaN, and MATLAB/Python/R/Julia behave differently in this
        ! respect. Moreover, MATLAB defines max(X) = [] if X == [], differing from mathematics
        ! and other languages.
        gamma = ZERO  ! The steplength of the projection step to be taken.
        if (any(resact(1:nact) > 1.0E-4_RP * delta)) then
            ! Set DPROJ to the shortest move (projection step) from S to the boundaries of the
            ! active constraints. We will use DPROJ to modify PSD.
            dproj = matprod(qfac(:, 1:nact), solve(transpose(rfac(1:nact, 1:nact)), resact(1:nact)))
            !!MATLAB: dproj = qfac(:, 1:nact) * (rfac(1:nact, 1:nact)' \ resact(1:nact))

            ! The vector DPROJ is also the shortest move from S + PSD to the boundaries of the
            ! active constraints (this is because PSD is parallel to the boundaries of the active
            ! constraints). Set GAMMA to the greatest steplength of this move that satisfies both
            ! the trust region bound and the linear constraints.
            ds = inprod(dproj, s + psd)
            dd = sum(dproj**2)
            resid = delsq - sum((s + psd)**2)
            ! Powell's condition for the following IF: RESID > 0.
            if (resid > 0 .and. dd > EPS * delsq .and. .not. is_nan(ds)) then
                ! Set GAMMA to the greatest value so that S + PSD + GAMMA*DPROJ satisfies the trust
                ! region bound. SQRTD: square root of a discriminant. Powell's code for SQRTD is
                ! SQRT(DS * DS + DD * RESID), which may be below ABS(DS) due to underflow in DS*DS.
                sqrtd = maxval([sqrt(ds * ds + dd * resid), abs(ds), sqrt(dd * resid)])
                if (ds <= 0) then
                    gamma = (sqrtd - ds) / dd
                else
                    gamma = resid / (sqrtd + ds)
                end if
                ! GAMMA < 0 should not happen. GAMMA can be 0 or NaN when, e.g., DS or DD becomes
                ! Inf. Powell's code does not handle this.
                if (gamma < 0 .or. .not. is_finite(gamma)) then
                    gamma = 0
                end if

                ! Reduce GAMMA so that the move along DPROJ also satisfies the linear constraints.
                ad = -ONE
                ad(trueloc(resnew > 0)) = matprod(dproj, amat(:, trueloc(resnew > 0)))
                frac = ONE
                restmp(trueloc(ad > 0)) = resnew(trueloc(ad > 0)) - matprod(psd, amat(:, trueloc(ad > 0)))
                frac(trueloc(ad > 0)) = restmp(trueloc(ad > 0)) / ad(trueloc(ad > 0))
                gamma = minval([gamma, ONE, frac])  ! GAMMA = MINVAL([GAMMA, ONE, FRAC(TRUELOC(AD>0))])
            end if
        end if

        ! Set the next direction for seeking a reduction in the model function subject to the trust
        ! region bound and the linear constraints.
        ! Do NOT write D = PSD + GAMMA*DPROJ, as DPROJ may contain NaN/Inf, in which case GAMMA = 0.
        if (gamma > 0) then
            d = psd + gamma * dproj  ! Modified searching direction.
            itercg = -1
        else
            d = psd  ! Original searching direction.
            itercg = 0
        end if
    end if
    itercg = itercg + 1_IK
    ! After the above line, ITERCG = 0 iff GETACT has been just called, and D is not PSD but a
    ! modified step.

    ! Set ALPHA to the steplength from S along D to the trust region boundary. Return if the first
    ! derivative term of this step is sufficiently small or if no further progress is possible.
    resid = delsq - ss
    dg = inprod(d, g)
    ds = inprod(d, s)
    dd = inprod(d, d)
    ! Powell's condition for the following IF: (RESID <= 0 .OR. DG >= 0). If DD is tiny (so is DS),
    ! ALPHA may be mistakenly calculated as a huge value due to rounding errors, as observed on
    ! 20221205. Therefore, we exit when DD is small. The test for DG is covered by the IF after the
    ! calculation of ALPHA.
    if (resid <= 0 .or. dd <= EPS * delsq .or. is_nan(ds)) then
        exit
    end if
    ! SQRTD: square root of a discriminant. Powell's code for SQRTD is SQRT(DS * DS + DD * RESID),
    ! which may be below ABS(DS) due to underflow in DS*DS.
    sqrtd = maxval([sqrt(ds * ds + dd * resid), abs(ds), sqrt(dd * resid)])
    if (ds <= 0) then
        alpha = (sqrtd - ds) / dd
    else
        alpha = resid / (sqrtd + ds)
    end if
    ! ALPHA < 0 should not happen. ALPHA can be 0 or NaN when, e.g., DS or DD becomes Inf. Powell's
    ! code does not handle this.
    if (alpha <= 0 .or. .not. is_finite(alpha)) then
        exit
    end if

    ! Powell's condition for the following IF: -ALPHA * DG <= TOL * REDUCT. Note that the EXIT
    ! will be triggered if DG >= 0, as ALPHA >= 0.
    if (-alpha * dg <= tol * reduct .or. is_nan(alpha * dg)) then
        exit
    end if

    ! Set DHD to the curvature of the model along D. Then reduce ALPHA if necessary to the value
    ! that minimizes the model.
    hd = matprod(hess, d)
    dhd = inprod(d, hd)
    alpht = alpha
    if (dg + alpha * dhd > 0) then
        alpha = -dg / dhd
    end if

    ! Make a further reduction in ALPHA if necessary to preserve feasibility.
    alphm = alpha
    ad = -ONE
    ad(trueloc(resnew > 0)) = matprod(d, amat(:, trueloc(resnew > 0)))
    frac = alpha
    frac(trueloc(ad > 0)) = resnew(trueloc(ad > 0)) / ad(trueloc(ad > 0))
    frac(trueloc(is_nan(frac))) = alpha
    jsav = 0
    if (any(frac < alpha)) then
        jsav = int(minloc(frac, dim=1), kind(jsav))
        alpha = frac(jsav)
    end if
    !----------------------------------------------------------------------------------------------!
    ! Alternatively, JSAV and ALPHA can be calculated as below.
    ! !JSAV = INT(MINLOC([ALPHA, FRAC], DIM=1), KIND(JSAV)) - 1_IK
    ! !ALPHA = MINVAL([ALPHA, FRAC])  ! This line cannot be exchanged with the last.
    ! We prefer our implementation as the code is more explicit; in addition, it is more flexible:
    ! we can change the condition ANY(FRAC < ALPHA) to ANY(FRAC < (1 - EPS) * ALPHA) or
    ! ANY(FRAC < (1 + EPS) * ALPHA), depending on whether we believe a false positive or a false
    ! negative of JSAV > 0 is more harmful.
    !----------------------------------------------------------------------------------------------!

    ! Post-process ALPHA according to some prior information.
    ! N.B.:
    ! 1. Since we set ALPHA=1 when ITERCG=0, the ALPHA calculated above is needed only if ITERCG>0.
    ! 2. Zaikun 20220821: In theory, shouldn't this post-processing change nothing? According to
    ! a test on 20220821, it does change ALPHA sometimes. Strange! Why?
    if (itercg == 0) then  ! Iff GETACT has been called, and D is not PSD but a modified step.
        ! By the definition of D, ALPHA = ONE is the largest ALPHA so that S + ALPHA*D satisfies the
        ! linear and trust region constraints.
        alpha = ONE
    elseif (itercg == 1 .and. gamma <= 0) then  ! Iff GETACT has been called, and D is not modified.
        ! Due to the scaling of PSD, S + D satisfies the linear and trust region constraints.
        alpha = max(alpha, ONE)
    else
        alpha = max(alpha, ZERO)
    end if

    ! Set ALPHA to the minimum between ALPHA and ALPHM, namely the steplength obtained by minimizing
    ! the quadratic model along D.
    alpha = min(alpha, alphm)

    ! Update S, G.
    sold = s
    s = s + alpha * d
    ss = sum(s**2)
    if (.not. is_finite(ss)) then
        s = sold
        exit
    end if
    g = g + alpha * hd
    if (.not. is_finite(sum(abs(g)))) then
        exit
    end if

    ! Update RESNEW.
    restmp = resnew - alpha * ad  ! Only RESTMP(TRUELOC(RESNEW > 0)) is needed.
    resnew(trueloc(resnew > 0)) = max(TINYCV, restmp(trueloc(resnew > 0)))
    !!MATLAB: mask = (resnew > 0); resnew(mask) = max(TINYCV, resnew(mask) - alpha * ad(mask));

    ! Update RESACT. This is done iff GETACT has been called, and D is not PSD but a modified step.
    !----------------------------------------------------------------------------------------------!
    ! Zaikun 20220821: There seems be a typo here. Powell's original code does not take ALPHA into
    ! account. Then RESACT seems to correspond to S + D, where D is defined as PSD + GAMMA*DPROJ
    ! during the modification procedure after GETACT is called. Without this modification, RESACT
    ! would remain unchanged because D = PSD, which is in the null space of the active constraints.
    ! The GAMMA*DPROJ component in the modified step D reduces RESACT by GAMMA*RESACT. However,
    ! since S is updated to S + ALPHA*D, shouldn't RESACT be reduced by ALPHA*GAMMA*RESACT?
    ! Note that Powell chose to update RESACT after ALPHA is calculated (instead of right after
    ! GAMMA is calculated), which might be an indication that he wanted to take ALPHA into account.
    ! In the following code, we try correcting this apparent typo, but it has little impact on the
    ! performance of LINCOA according to a test on 20220821.
    if (itercg == 0) then
        resact(1:nact) = (ONE - alpha * gamma) * resact(1:nact)
        !resact(1:nact) = (ONE - gamma) * resact(1:nact)  ! Powell's code.
    end if
    !----------------------------------------------------------------------------------------------!

    ! Update REDUCT, the reduction up to now.
    reduct = reduct - alpha * (dg + HALF * alpha * dhd)
    if (reduct <= 0 .or. is_nan(reduct)) then
        s = sold
        exit
    end if

    ! Test for termination.
    if (alpha >= alpht .or. -alphm * (dg + HALF * alphm * dhd) <= tol * reduct) then
        exit
    end if

    ! Branch to a new loop if there is a new active constraint.
    ! When JSAV > 0, Powell's code branches back with NEWACT = .TRUE. only if ||S|| <= 0.8*DELTA,
    ! and it exits if ||S|| > 0.8*DELTA, as mentioned at the end of Section 3 of Powell 2015. The
    ! motivation seems to avoid small steps that changes the active set, because GETACT is expensive
    ! in flops. However, according to a test on 20220820, removing this condition (essentially
    ! replacing it with ||S|| < DELTA) improves the performance of LINCOA a bit. This may lead to
    ! small steps, but tiny steps will lead to tiny reductions and trigger an exit.
    newact = (jsav > 0)
    if (newact) then
        cycle
    end if

    ! If N-NACT CG iterations has been taken in the current null space (corresponding to the
    ! current "active set"), then, in theory, a stationary point in this subspace has been found.
    ! If the "active set" is the true active set, then a stationary point of the
    ! linearly-constrained trust region subproblem is found. So a termination is reasonable.
    ! However, the "active set" is not precisely the true active set, is it? See (3.5) of Powell
    ! (2015) and the comments at the beginning of the GETACT subroutine. Also, we should take into
    ! account the modification after GETACT is called.
    if (itercg >= n - nact) then  ! ITERCG > N - NACT is impossible.
        exit
    end if

    ! Calculate the next search direction, which is conjugate to the previous one if ITERCG /= NACT.
    ! N.B.: NACT < 0 is impossible unless GETACT is buggy; NACT = 0 can happen, particularly if
    ! there is no constraint. In theory, the code for the second case below covers the first as well.
    if (nact <= 0) then
        pg = g
    else
        pg = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
        !!MATLAB: pg = qfac(:, nact+1:n) * (g' * qfac(:, nact+1:n))';
    end if

    if (itercg == 0) then  ! Iff GETACT has been called, and D is not PSD but a modified step.
        beta = ZERO
    else
        beta = inprod(pg, hd) / dhd
    end if
    d = -pg + beta * d
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(s) == n .and. all(is_finite(s)), 'SIZE(S) == N, S is finite', srname)
    ! Due to rounding, it may happen that ||S|| > DELTA, but ||S|| > 2*DELTA is highly improbable.
    call assert(norm(s) <= TWO * delta, '||S|| <= 2*DELTA', srname)
end if

end subroutine trslin
!--------------------------------------------------------------------------------------------------!
! Zaikun 20220417:
! For PG, the schemes below work evidently worse than the one above in a test on 20220417. Why?
!-----------------------------------------------------------------------!
! VERSION 1:
! !pg = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
!-----------------------------------------------------------------------!
! VERSION 2:
! !if (2 * nact < n) then
! !    pg = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
! !else
! !    pg = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
! !end if
!-----------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


end submodule trslin_mod