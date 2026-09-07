submodule (trustregion_mod) trsunc_mod
!--------------------------------------------------------------------------------------------------!
! Global subproblem solver for unconstrained trust-region subproblem:
! 
! 	min_{x} Q(x) := dot(g, x) + 0.5 * x^T * H * x, subject to ||x|| <= delta
! 
! The solver can also handle the adaptive cubic regularization subproblem,
! 	min_{x} Q(x) + (1/3) * delta * ||x||^3 = dot(g, x) + 0.5 * x^T * H * x + (1/3) * delta * ||x||^3
! which is globally optimized using very similar methods.
! 
! More details about these problems are available in:
!   [CGT2000] A. R. Conn, N. I. M. Gould, P. L. Toint. Trust-Region Methods. SIAM (2000).
!   [CGT2022] C. Cartis, N. I. M. Gould, P. L. Toint. Evaluation Complexity of Algorithms for Nonconvex Optimization. SIAM (2022).
! 
! The algorithm implemented here is primarily based on:
!   [GRT2010] N. I. M. Gould, D. P. Robinson, H. S. Thorne. On solving trust-region and other regularised subproblems in optimization.
!             Mathematical Programming Computation 2:1 (2010), pp. 21-57. 
! The same algorithm is also implemented in routines TRS/RQS in GALAHAD; some parameter values are taken from that code.
! The key ideas of [GRT2010] are presented in summarized form in Chapter 9 of [CGT2022].
! 
! This implementation has a fallback strategy: if the main algorithm fails, use the Cauchy step to ensure
! sufficient decrease is achieved. This is a step of the form x = -alphaC*g, for some alphaC>0 corresponding to
! an exact linesearch.
!--------------------------------------------------------------------------------------------------!

implicit none

contains

subroutine trglob(delta, g_in, hess_in, lambda, x, is_tr)
	! Common modules
	use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, HALF, REALMIN, ZERO, TENTH, EPS, DEBUGGING
	use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_inf
	use, non_intrinsic :: linalg_mod, only : inprod, issymmetric, norm, project, matprod
	use, non_intrinsic :: univar_mod, only : circle_min

	use, non_intrinsic :: qdec_mod, only : quadform, qdec, crdec

	implicit none

	! Inputs
	real(RP), intent(in) :: delta
	real(RP), intent(in) :: g_in(:)  ! G_IN(N)
	real(RP), intent(in) :: hess_in(:, :)  ! HESS_IN(N, N)
	logical, intent(in) :: is_tr

	! Outputs
	real(RP), intent(out) :: lambda
	real(RP), intent(out) :: x(:)  ! X(N)

	! Locals
	character(len=*), parameter :: srname = 'TRGLOB'
	integer(IK) :: n
	integer(IK) :: status
	logical :: scaled
	logical :: result
	real(RP) :: g_scal(size(g_in))
	real(RP) :: H_scal(size(hess_in, 1), size(hess_in, 2))
	real(RP) :: H_plus_lambda_I(size(hess_in, 1), size(hess_in, 2))
	real(RP) :: modscal
	real(RP) :: lambdaL, lambdaU, lambdaC
	real(RP) :: dval ! output from safe_cholesky
	real(RP) :: norm_check  ! desired value of ||x(lambda)||, potentially updated at each iteration
	real(RP) :: rho
	real(RP) :: gnorm
	real(RP) :: cos_z_g
	logical :: found_lambdaC_in_L
	logical :: potential_hard_case
	logical :: hard_case
	integer(IK) :: it
	integer(IK) :: i
	integer(IK) :: it1, it2
	real(RP) :: lambda1_neg1, lambda2_2, lambda3_2, lambdaT
	real(RP) :: small_width
	integer(IK) :: nk
	real(RP) :: gammak, new_lambda
	real(RP) :: cauchy_decrease, current_decrease
	real(RP) :: alpha
	logical :: use_cauchy_step
	real(RP) :: xcauchy(size(g_in))
	real(RP) :: z(size(g_in))
	real(RP) :: lambda_plus, lambda1_neg1, lambda3_2
	real(RP) :: old_lambdaC
	real(RP) :: lambda_width
	logical, parameter :: verbose = .true.

	! Solver parameters
	real(RP), parameter :: gamma = 1.0  ! for setting initial lambda (eq 3.51), value taken from GALAHAD/trs.f90
	real(RP), parameter :: theta = 0.01  ! for setting initial lambda (eq 3.51), value taken from GALAHAD/trs.f90
	real(RP), parameter :: z_parallel_g_thresh = 1e-8  ! identify when z is parallel to g (for inverse iterations)
	real(RP), parameter :: bisection_thresh = 1e-12  ! identify when to terminate main bisection loop with failure
	real(RP), parameter :: small_width_thresh = 1e-5  ! if lambdaC in G, update to lambdaT if not too close to the boundary
	integer(IK), parameter :: num_inverse_iters_lambda_in_G = 1  ! number of inverse iterations to do when lambdaC in G
	real(RP), parameter :: omega = 1.0  ! not sure what a good value is for this (Algorithm 3.3 updating ratio)
	integer(IK), parameter :: num_potential_hard_case_iters = 10  ! number of iterations of Algorithm 3.3
	real(RP), parameter :: lambda_convergence_thresh = 1e-12  ! how close successive iterates of Algorithm 3.3 are before termination
	real(RP), parameter :: potential_hard_case_thresh_decrease = 0.1  ! how much to decrease potential_hard_case_thresh when it's too large
	integer(IK), parameter :: num_refinement_iters = 20  ! max number of refinement iterations
	real(RP), parameter :: boundary_thresh = 1e-12  ! how close should ||x|| be to delta before terminating
	real(RP), parameter :: rounding_perturbation = 1e-12  ! how much to perturb final lambdaC by to ensure it's just above -lambda1
	real(RP) :: potential_hard_case_thresh = 1e-2  ! interval width when to check for hard case (will be reduced if necessary)

	! Sizes
	n = int(size(g_in), kind(n))

	! Preconditions
	if (DEBUGGING) then
		call assert(n >= 1, 'N >= 1', srname)
		call assert(delta > 0, 'DELTA > 0', srname)
		call assert(size(g_in) == n, 'SIZE(G) = N', srname)
		call assert(size(hess_in, 1) == n .and. issymmetric(hess_in), 'HESS is an NxN symmetric matrix', srname)
		call assert(size(x) == n, 'SIZE(X) == N', srname)
	end if

	if (is_tr .and. delta < boundary_thresh) then 
		! Easy quit: if constraint is ||s|| <= 0, then only feasible solution is s=0
		x = ZERO
		lambda = ZERO
		return
	end if

	if (maxval(abs(g_in)) > 1.0E12) then   ! The threshold is empirical.
		modscal = max(TWO * REALMIN, ONE / maxval(abs(g_in)))  ! MAX: precaution against underflow.
		g_scal = g_in * modscal
		H_scal = hess_in * modscal
		scaled = .true.
	else
		modscal = ONE  ! This value is not used, but Fortran compilers may complain without it.
		g_scal = g_in
		H_scal = hess_in
		scaled = .false.
	end if

	result = .false.  ! have we found a value for x yet?

	! Find initial lambda region and current guess
	status = 0

	call initlda(H_scal, g_scal, delta, lambdaL, lambdaU, is_tr)
	if (lambdaL == ZERO) then
		lambdaC = ZERO
	else
		lambdaC = max(gamma * sqrt(lambdaL * lambdaU), lambdaL + theta * (lambdaU - lambdaL))
	end if

	if (verbose) then
		print *, "Initially, lambdaL = ", lambdaL, " and lambdaU = ", lambdaU
	end if

	! ---------------- Main loop ---------------- !

	! For (TRS), check for interior solution before proceeding further
	if (is_tr .and. lambdaL == ZERO) then
		call solvekkt(H_scal, g_scal, ZERO, x, dval, H_plus_lambda_I, status)
		if (status == 0 .and. sqrt(sum(x**2)) <= delta + boundary_thresh * max(ONE, delta)) then
			if (verbose) then
				print *, "Interior solution found, terminating"
			end if
			return
		end if
	end if

	! Initialize z to a unit vector orthogonal to g (used for all inverse iterations)
	z = ONE / sqrt(real(n, RP))

	! First, make sure z is not parallel to g...
	gnorm = sqrt(sum(g_scal**2))
	if (gnorm .ne. ZERO) then
		! If g=0, then z is definitely not parallel to g
		cos_z_g = inprod(z, g_scal) / gnorm
		if (abs(cos_z_g) >= ONE - z_parallel_g_thresh) then
			! cos(z,g) ~ 1, so g is basically a vector with all entries the same...
			z(1) = -z(1)  ! this will definitely make z not parallel to g
		end if
		! Now make z orthogonal to g using the formula z = z - (ghat.T * z) * ghat, where ghat = g / ||g||
		z(1:n) = z(1:n) - inprod(z, g_scal) / sum(g_scal**2) * g_scal  ! z += [-dot(g,z)/||g||^2] * g
	end if

	! *** Phase 1 - find a good estimate of lambdaC ***
	!
	! Regions of lambda for (TRS):
	!   - N = H(lambda) is negative semidefinite --> lambda in (-inf, lambdaS]
	!   - L = H(lambda) is positive definite and ||x(lambda)|| >= delta, --> lambda in (lambdaS, lambda*]
	!   - G = H(lambda) is positive definite and ||x(lambda)|| < delta --> lambda in (lambda*, inf)
	!
	! For (ARC), same except check for L vs G is if ||x(lambda)|| >= lambda/delta
	! As lambda -> lambdaS from above (i.e. from L), ||x(lambda)|| -> +infty
	! At lambda = lambda* (i.e. on L/G boundary), we have ||x(lambda)|| = delta or lambda/delta
	! As lambda -> +infty (i.e. in G), ||x(lambda)|| -> 0
	!
	! where H(lambda)=H+lambda*I and x(lambda) solves H(lambda)*x(lambda)=-g

	found_lambdaC_in_L = .false.
	potential_hard_case = .false.
	hard_case = .false.
	it = 0

	lambda_width = max(bisection_thresh * max(abs(lambdaL), abs(lambdaU)), bisection_thresh)
	do while (abs(lambdaU - lambdaL) > lambda_width .and. (.not. found_lambdaC_in_L) .and. (.not. hard_case))
		if (lambdaL >= lambdaU) then
			exit  ! always expect lambdaL < lambdaU, so stop the loop if this is violated
		end if

		if (verbose) then
			print *, "It", it, ", lambdaL = ", lambdaL, "lambdaC = ", lambdaC, ", lambdaU = ", lambdaU
		end if

		it = it + 1

		if (.not. potential_hard_case) then
			! Identify which region lambdaC is in
			call solvekkt(H_scal, g_scal, lambdaC, x, dval, H_plus_lambda_I, status)

			if (status == 0) then
				! H + lambdaC*I is positive definite
				if (is_tr) then
					norm_check = delta
				else
					norm_check = lambdaC / delta
				end if
				if (sqrt(sum(x**2)) >= norm_check) then
					! lambdaC in L, i.e. ||x(lambda)|| >= delta or lambdaC/delta
					! This is the good region, where we switch to the fast-converging phase 2
					found_lambdaC_in_L = .true.
					if (verbose) then
						print *, "lambdaC in L, stopping"
					end if
					exit
				else
					! lambda in G, i.e. ||x(lambda)|| < delta (lambda too large)
					lambdaU = min(lambdaU, lambdaC)

					! Update lambdaL using inverse iteration
					do i = 1, num_inverse_iters_lambda_in_G
						! H_plus_lambda_I has already set and factorized in solvekkt above
						call cholsolve(H_plus_lambda_I, z)  ! z <-- (H+lambdaC*I) \ z
						z(1:n) = z(1:n) / sqrt(sum(z**2))  ! z <-- z / ||z||
					end do

					rho = quadform(H_scal, z)
					lambdaL = max(lambdaL, -rho)

					! Get under-approximations of lambdaC using Taylor approximations, and update lambdaL based on this
					lambdaT = lambdaL
					call newlda(H_plus_lambda_I, x, lambdaC, delta, -ONE, 1, lambda1_neg1, is_tr, status)
					if (status == 0) then
						lambdaT = max(lambdaT, lambda1_neg1)
					end if
					call newlda(H_plus_lambda_I, x, lambdaC, delta, TWO, 2, lambda2_2, is_tr, status)
					if (status == 0) then
						lambdaT = max(lambdaT, lambda2_2)
					end if
					call newlda(H_plus_lambda_I, x, lambdaC, delta, TWO, 3, lambda3_2, is_tr, status)
					if (status == 0) then 
						lambdaT = max(lambdaT, lambda3_2)
					end if

					small_width = small_width_thresh * abs(lambdaU - lambdaL)
					if ((lambdaT >= lambdaL + small_width) .and. (lambdaT <= lambdaU - small_width)) then
						lambdaC = lambdaT
					else
						lambdaC = max(gamma * sqrt(lambdaL * lambdaU), lambdaL + theta * (lambdaU - lambdaL))
					end if
				end if
			else if (status > 0) then
				! lambda in N, i.e. H+lambdaC*I is not positive definite (lambdaC too small)
				! Here, x and d are set so that x.T * (H + lambdaC*I + d * ek * ek.T) * x = 0
				lambdaL = max(lambdaL, lambdaC)
				lambdaL = max(lambdaL, dval / sum(x**2) + lambdaC)
				lambdaC = max(gamma * sqrt(lambdaL * lambdaU), lambdaL + theta * (lambdaU - lambdaL))
			else
				! Error in solvekkt
				if (verbose) then
					print *, "Error in solvekkt, stopping"
				end if
				exit
			end if

			! Check if we are in the potential hard case now
			if (abs(lambdaU - lambdaL) <= max(potential_hard_case_thresh * max(abs(lambdaL), abs(lambdaU)), potential_hard_case_thresh)) then
				! Identified potential hard case
				if (verbose) then
					print *, "Identified potential hard case, lambdaL = ", lambdaL, ", lambdaU = ", lambdaU
				end if
				potential_hard_case = .true.
			else 
				potential_hard_case = .false.
			end if
			! end regular bisection step
		else
			! Potential hard case
			if (verbose) then
				print *, "Potential hard case"
			end if
			old_lambdaC = lambdaC
			lambdaC = lambdaU  ! definitely an overestimate of -lambda1

			do it2 = 1, num_potential_hard_case_iters
				if (verbose) then
					print *, "It2 = ", it2, ", trying lambdaC = ", lambdaC
				end if
				call solvekkt(H_scal, g_scal, lambdaC, x, dval, H_plus_lambda_I, status)

				if (status == 0) then
					! H + lambdaC*I is positive definite (expected)
					if (is_tr) then
						norm_check = delta
					else
						norm_check = lambdaC / delta
					end if
					if (sqrt(sum(x**2)) >= norm_check) then
						! lambdaC in L, i.e. ||x(lambda)|| >= delta or lambdaC/delta
						! This is the good region, where we switch to the fast-converging phase 2
						found_lambdaC_in_L = .true.
						if (verbose) then
							print *, "lambdaC in L, stopping"
						end if
						exit
					end if
					! Otherwise, update lambdaC using inverse iteration
					if (it > 5) then
						nk = 1  ! don't need too many iterations
					else 
						nk = 2
					end if
					do i = 1, nk
						! H_plus_lambda_I has already set and factorized in solvekkt above
						call cholsolve(H_plus_lambda_I, z)  ! z <-- (H+lambdaC*I) \ z
						z(1:n) = z(1:n) / sqrt(sum(z**2))  ! z <-- z / ||z||
					end do

					! Update lambdaC
					rho = quadform(H_scal, z)
					if (nk == 1) then
						gammak = 1.5
					else
						gammak = 3.0
					end if
					new_lambda = -rho + omega * (lambdaC + rho) ** (gammak)
					if (abs(new_lambda - lambdaC) < lambda_convergence_thresh) then
						hard_case = .true.
						exit
					end if
					if (new_lambda < lambdaC) then
						lambdaC = new_lambda
					else
						! If we start lambdaC sufficiently close to a good value, this should never happen
						! So, stop the 'potential hard case' iteration and go back to regular bisection phase
						if (verbose) then
							print *, "Started nearly hard case too soon, trying again"
						end if
						potential_hard_case = .false.
						lambdaC = old_lambdaC
						potential_hard_case_thresh = potential_hard_case_thresh * potential_hard_case_thresh_decrease
						exit
					end if
				else if (status > 0) then
					! H + lambdaC*I is not positive definite
					! This happens when lambdaC is slightly smaller than -lambda_{max}(H), due to rounding errors
					if (verbose) then
						print *, "Potential hard case failure (lambdaC too small from rounding errors) - treat as hard case"
					end if
					hard_case = .true.
					exit
				else
					! Error in solvekkt, stopping
					if (verbose) then
						print *, "Error in solvekkt, stopping"
					end if
					exit
				end if
			end do  ! end potential hard case iteration loop

			! Stop the bisection loop
			if (found_lambdaC_in_L .or. hard_case .or. status < 0) then
				exit
			end if

			! If the potential hard case didn't find lambdaC in L, then lambdaC is still in G and is a potential upper bound
			lambdaU = min(lambdaC, lambdaU)
		end if ! end potential hard case

		! Update desired interval width for termination at start of next iteration
		lambda_width = max(bisection_thresh * max(abs(lambdaL), abs(lambdaU)), bisection_thresh)
	end do ! end bisection loop

	! End of main lambdaC search phase
	! Three possibilities to get here:
	! 1. Found lambdaC in L, which needs to be increased via a refinement process
	! 2. In the hard case, which has a solution with a special form
	! 3. Error in above search, terminate with failure

	if (found_lambdaC_in_L) then
		! Refine an estimate lambdaC in L
		if (verbose) then
			print *, "Phase 2: refining lambdaC = ", lambdaC
		end if
	
		do it2 = 1, num_refinement_iters
			call newlda(H_plus_lambda_I, x, lambdaC, delta, -ONE, 1, lambda1_neg1, is_tr, status);
			if (status == 0) then
				call newlda(H_plus_lambda_I, x, lambdaC, delta, TWO, 3, lambda3_2, is_tr, status);
				if (status == 0) then
					lambda_plus = max(lambda1_neg1, lambda3_2)
				else
					if (verbose) then
						print *, "newlda2 failed"
					end if
					lambda_plus = lambda1_neg1
				end if
			else
				if (verbose) then
					print *, "newlda1 failed"
				end if
				lambda_plus = lambdaC
			end if

			if (verbose) then
				print *, "Refining iteration", it2, " found new lambdaC = ", lambda_plus
			end if
			
			if (abs(lambdaC - lambda_plus) < EPS * max(1.0, abs(lambdaC))) then
				! termination from GALAHAD/trs.f90 -- refinement iteration not achieving much
				lambdaC = lambda_plus
				if (verbose) then
					print *, "Terminating refinement phase (limited progress)"
				end if
				exit
			end if
			lambdaC = lambda_plus

			! Recompute factorization
			call solvekkt(H_scal, g_scal, lambdaC, x, dval, H_plus_lambda_I, status)
			if (status == 0) then
				! H + lambdaC*I is positive definite
				! Check if we are near the desired norm
				if (is_tr) then
					norm_check = delta
				else
					norm_check = lambdaC / delta
				end if
				if (abs(sqrt(sum(x**2)) - norm_check) < boundary_thresh * max(ONE, norm_check)) then
					! Terminating near boundary (success)
					if (verbose) then
						print *, "Terminating near boundary, success"
					end if
					exit
				end if
			else
				! In the refinement phase (lambdaC in L), we should never be able to produce a new lambdaC in N
				! i.e. H + lambdaC*I should always be positive definite

				! status > 0 --> H + lambdaC*I is not positive definite (left the good region)
				if (verbose) then
					if (status > 0) then
						print *, "lambdaC has left the good region!"
					else
						print *, "Positive definite Cholesky solve failed, stopping"
					end if
				end if
				! status <= 0 --> error in solvekkt
				exit
			end if
		end do

		! End of refinement phase, terminate with good estimate...
		if (status == 0) then
			! Final solve with lambdaC
			if (verbose) then
				print *, "Final solve with lambdaC = ", lambdaC
			end if
			call solvekkt(H_scal, g_scal, lambdaC, x, dval, H_plus_lambda_I, status)
			if (status == 0) then
				if (is_tr) then
					norm_check = delta
				else
					norm_check = lambdaC / delta
				end if
				if (abs(sqrt(sum(x**2)) - norm_check) < boundary_thresh * max(ONE, norm_check)) then
					result = .true.
					if (verbose) then
						print *, "Phase 2 success"
					end if
				else
					! Scale x to have the desired norm
					x(1:n) = x(1:n) * (norm_check / sqrt(sum(x**2)))
					result = .true.
					if (verbose) then
						print *, "Phase 2 solution too far from boundary, rescaling"
					end if
				end if
			else
				! status > 0 --> Easy case error: final lambdaC gave indefinite Hessian
				! status <= 0 --> Positive definite Cholesky solve failed
				if (verbose) then
					if (status > 0) then
						print *, "Easy case error: final lambdaC gave indefinite Hessian"
					else 
						print *, "Positive definite Cholesky solve failed, stopping"
					end if
				end if
				result = .false.
			end if
		else
			! Error in refinement phase
			if (verbose) then
				print *, "Error in refinement phase"
			end if
			result = .false.
		end if
	else if (hard_case) then
		! Hard case
		if (verbose) then
			print *, "Hard case"
		end if
		! Here, z is a good estimate of u1, a unit eigenvector corresponding to lambda1
		lambdaC = -quadform(H_scal, z)  ! lambdaC = -lambda1
		call solvekkt(H_scal, g_scal, lambdaC, x, dval, H_plus_lambda_I, status)

		if (status > 0) then
			! H + lambdaC*I not positive definite, try again with slightly larger lambdaC to avoid rounding errors
			lambdaC = lambdaC + max(rounding_perturbation, rounding_perturbation * abs(lambdaC))
			call solvekkt(H_scal, g_scal, lambdaC, x, dval, H_plus_lambda_I, status)
		end if

		if (status == 0) then
			! Final solution is x + alpha*z, with alpha chosen to give correct vector norm
			call hardstep(x, z, delta, lambdaC, alpha, is_tr, status)
			if (status == 0) then
				if (verbose) then
					print *, "Using alpha = ", alpha
				end if
				x(1:n) = x(1:n) + alpha * z(1:n)
			end if
			! status != 0, i.e. if no roots to the quadratic, then ||xs|| sufficiently large already, so nothing to do
			result = .true.
		else
			! status > 0 --> Rounding errors, Rayleigh quotient gave eigenvalue underestimate
			! status <= 0 --> Positive definite Cholesky solve failed
			if (verbose) then
				if (status > 0) then
					print *, "Rounding errors: Rayleigh quotiend gave eigenvalue underestimate"
				else
					print *, "Positive definite Cholesky solve failed"
				end if
			end if
			result = .false.
		end if
		! end hard case
	else
		! Error in bisection phase
		if (verbose) then
			print *, "Error in bisection phase"
		end if
		result = .false.
	end if

	! Catch any inf/NaN issues here
	do i = 1, n 
		if (is_nan(x(i)) .or. is_inf(x(i))) then
			result = .false.
			exit
		end if
	end do

	! Compute Cauchy step as safeguard - use if above computation failed, or didn't get sufficient decrease
	call cauchy(H_scal, g_scal, delta, is_tr, xcauchy)  ! set xcauchy to Cauchy step
	use_cauchy_step = .false.

	if (result) then
		! If global minimizer computation above succeeded, make sure got at least Cauchy decrease
		use_cauchy_step = .false.
		if (is_tr) then
			cauchy_decrease = qdec(g_scal, H_scal, xcauchy)
			current_decrease = qdec(g_scal, H_scal, x)
		else
			cauchy_decrease = crdec(g_scal, H_scal, xcauchy, delta)
			current_decrease = crdec(g_scal, H_scal, x, delta)
		end if
		
		if (cauchy_decrease > current_decrease) then
			! Global step didn't achieve sufficient decrease, using Cauchy step instead
			if (verbose) then
				print *, "Global step didn't achieve sufficient decrease, using Cauchy step instead"
			end if
			use_cauchy_step = .true.
		else if (is_tr .and. (sqrt(sum(x**2)) > delta + boundary_thresh * max(1.0, delta))) then
			! Global step outside feasible region, using Cauchy step instead
			if (verbose) then
				print *, "Global step outside feasible region, using Cauchy step instead"
			end if
			use_cauchy_step = .true.
		else
			use_cauchy_step = .false.
		end if
	else
		! Global step calculation failed, using Cauchy step instead
		if (verbose) then
			print *, "Global step calculation failed, using Cauchy step instead"
		end if
		use_cauchy_step = .true.
	end if

	! Set final lambda value
	lambda = lambdaC

	if (use_cauchy_step) then 
		x(1:n) = xcauchy(1:n)
		lambda = -ONE  ! flag failure
	end if

end subroutine trglob


subroutine initlda(H, g, delta, lambdaL, lambdaU, is_tr)
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : issymmetric

	implicit none

	! Inputs
	real(RP), intent(in) :: H(:,:)  ! H(N,N)
	real(RP), intent(in) :: g(:)  ! G(N)
	real(RP), intent(in) :: delta
	logical, intent(in) :: is_tr

	! Outputs
	real(RP), intent(out) :: lambdaL
	real(RP), intent(out) :: lambdaU
	
	! Locals
	character(len=*), parameter :: srname = 'INITLDA'
	integer(IK) :: i
	integer(IK) :: n
	integer(IK) :: nroots
	! expand [lambdaL, lambdaU] interval slightly, in case correct lambda is at endpoint
	real(RP), parameter :: expand_interval_thresh = 1.0e-5
	! Bounds on min/max eigenvalues of H, satisfying:
	! -lambda_min(H) <= lambda_min_bound <-- NOTE SIGN FLIP
	! lambda_max(H) <= lambda_max_bound
	real(RP) :: lambda_min_bound, lambda_max_bound
	real(RP) :: H_normF, H_normInf, H_min_diag
	! Use Gershgorin discs to get next estimates on min/max eigenvalues
	! gershgorin_lower_bound <= lambda(H) <= gershgorin_upper_bound
	real(RP) :: gershgorin_lower_bound, gershgorin_upper_bound
	real(RP) :: normg
	real(RP) :: l1, l2
	real(RP) :: current_sum
	real(RP) :: center, radius

	! Sizes.
    n = int(size(g), kind(n))

    ! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(g) == n, 'SIZE(G) == N', srname)
        call assert(size(H, 1) == n .and. issymmetric(H), 'HESS is n-by-n and symmetric', srname)
    end if


	! Use Frobenius and infinity norms of H to get first estimates on min/max eigenvalues
	H_normF = sqrt(sum(H**2))

	H_normInf = sum(abs(H(1, :)))
	do i = 2, n  ! start from i=2
		current_sum = sum(abs(H(i, :)))
		if (current_sum > H_normInf) then
			H_normInf = current_sum
		end if
	end do

	H_min_diag = H(1, 1)
	do i = 2, n  ! start from i=2
		H_min_diag = min(H_min_diag, H(i,i))
	end do

	lambda_min_bound = min(H_normF, H_normInf)
	lambda_max_bound = min(H_normF, H_normInf)

	! Use Gershgorin discs to get next estimates on min/max eigenvalues
	! For each row, Gershgorin disc is B(diag value, sum abs off-diag values)
	center = H(1, 1)
	radius = sum(abs(H(1,:))) - abs(center)
	gershgorin_lower_bound = center - radius
	gershgorin_upper_bound = center + radius
	do i=2, n  ! start from i=2
		center = H(i, i)
		radius = sum(abs(H(i, :))) - abs(center)
		gershgorin_lower_bound = min(gershgorin_lower_bound, center - radius)
		gershgorin_upper_bound = max(gershgorin_upper_bound, center + radius)
	end do

	lambda_min_bound = min(lambda_min_bound, -gershgorin_lower_bound)  ! note sign flip
	lambda_max_bound = min(lambda_max_bound, gershgorin_upper_bound)

	! Norm of g
	normg = sqrt(sum(g**2))

	! Set bounds
	lambdaL = max(ZERO, -H_min_diag)
	lambdaU = ZERO
	if (is_tr) then
		lambdaL = max(lambdaL, normg / delta - lambda_max_bound)
		lambdaU = max(lambdaU, normg / delta + lambda_min_bound)
	else
		! Bounds given by largest root of: lambda^2 + lambda_{min/max} * lambda - ||g|| * delta = 0
		! Always has at exactly one positive root when ||g||*delta > 0
		! Always have delta>0 from earlier checks, so if ||g||=0 then either H convex -> lambdaC=0 -> x=0
		! or H indefinite -> hard case lambdaC=-lambda1
		! Recall lambda_min_bound = -lambda_min estimate
		if (normg == ZERO) then
			! lambdaL = -ve best estimate of lambda_min, or zero
			lambdaL = max(lambdaL, lambda_min_bound); 
			! lambdaU slightly above lambdaL
			lambdaU = lambdaL + max(expand_interval_thresh, expand_interval_thresh * lambdaL)
		else
			
			call quadroots(ONE, -lambda_min_bound, -normg * delta, l1, l2, nroots)
			if (nroots == 2) then
				lambdaU = max(lambdaU, max(l1, l2))
			end if

			call quadroots(ONE, lambda_max_bound, -normg * delta, l1, l2, nroots)
			if (nroots == 2) then
				lambdaL = max(lambdaL, max(l1, l2))
			end if
		end if
	end if

	! Just in case the true lambdaC is exactly equal to one of these, widen the interval slightly so it becomes an interior point
	if (lambdaL > ZERO) then
		lambdaL = lambdaL - max(expand_interval_thresh, expand_interval_thresh * lambdaL)
	end if
	! the above may make lambdaL slightly negative, which we don't want
	lambdaL = max(lambdaL, ZERO)
	if (lambdaU > ZERO) then 
		lambdaU = lambdaU + max(expand_interval_thresh, expand_interval_thresh * lambdaU)
	end if
end subroutine initlda

subroutine solvekkt(H, g, lambdaC, x, dval, H_plus_lambda_I, status)
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : issymmetric

	implicit none

	! Inputs
	real(RP), intent(in) :: H(:,:)  ! H(N,N)
	real(RP), intent(in) :: g(:)  ! G(N)
	real(RP), intent(in) :: lambdaC

	! Outputs
	real(RP), intent(out) :: x(:)  ! x(N)
	real(RP), intent(out) :: dval
	real(RP), intent(out) :: H_plus_lambda_I(:, :)  ! H_plus_lambda_I(N,N)
	integer(IK), intent(out) :: status

	! Locals
	character(len=*), parameter :: srname = 'SOLVEKKT'
	integer(IK) :: i
	integer(IK) :: n
	

	! Sizes.
    n = int(size(g), kind(n))

    ! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(g) == n, 'SIZE(G) == N', srname)
        call assert(size(H, 1) == n .and. issymmetric(H), 'HESS is n-by-n and symmetric', srname)
        call assert(size(x) == n, 'SIZE(X) == N', srname)
		call assert(size(H_plus_lambda_I, 1) == n .and. size(H_plus_lambda_I, 2) == n, 'H_p_l_I is n-by-n and symmetric', srname)
    end if

	H_plus_lambda_I = H 
	do i = 1, n 
		H_plus_lambda_I(i, i) = H_plus_lambda_I(i, i) + lambdaC
	end do

	call cholsafe(H_plus_lambda_I, x, dval, status)
	
	if (status == 0) then
		! H + lambdaC*I was positive definite: solve (H + lambdaC*I) * x = -g

		x = -ONE * g
		call cholsolve(H_plus_lambda_I, x)
	end if
end subroutine solvekkt

subroutine printmat(A)
	use, non_intrinsic :: consts_mod, only : RP, IK
    implicit none
    real(RP) :: A(:, :)
    integer :: i

    do i = 1, size(A, 1)
        ! "*(1x, i4)" means: repeat for all elements, add 1 space, use 4 spaces for integer
		print *, A(i, :)
    end do
end subroutine printmat

subroutine cholsafe(A, v, delta, status)
	! --------------------------------------------------------------- !
	! Safe in-place Cholesky factorization: A = L * L^T
	! 
	! Overwrites the upper and lower triangular parts of the matrix with L
	! (so that A remains symmetric, which makes cholsolve easier)
	! 
	! Outputs of this function are:
	! - v = vector of length n, updated if A is not positive definite
	! - delta = value updated if A is not positive definite
	! 
	! Return value 'status' is
	! - 0 = successful, A is positive definite [v, delta unchanged]
	! - Positive = leading minor of order (output value) is not 
	!              positive definite [v, delta changed]
	! 
	! If A is not positive definite (i.e. return value > 0), then v 
	! and delta are set so that v is a nonzero vector such that 
	!     v.T * (A + delta * ek * ek.T) * v = 0,
	! where k is the iteration when failure occurs.
	! --------------------------------------------------------------- !
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : issymmetric

	implicit none

	! Inputs/outputs
	real(RP), intent(inout) :: A(:,:)  ! A(N,N)
	real(RP), intent(out) :: v(:)  ! V(N)
	real(RP), intent(out) :: delta
	integer(IK), intent(out) :: status

	! Locals
	character(len=*), parameter :: srname = 'CHOLSAFE'
	integer(IK) :: i, j, k, l
	integer(IK) :: n
	real(RP) :: Lkk, sqrt_Lkk
	real(RP) :: adiag(size(A, 1))  ! adiag(N)

	! Sizes.
    n = int(size(A, 1), kind(n))

	! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(A, 1) == n .and. issymmetric(A), 'A is n-by-n and symmetric', srname)
		call assert(size(v) == n, 'SIZE(V) == N', srname)
    end if

	! Algorithm is from Section 7.3.7 of Conn, Gould, Toint, Trust-Region Methods, SIAM (2000)
	
	! Save diag(A) in VTMP, as we will need these values to compute delta if failure
	do i = 1, n 
		adiag(i) = A(i, i)
	end do

	! Set return values based on success, overridden if failure
	status = 0
	delta = ZERO
	v = ZERO

	do k = 1, n
		!print *, "cholsafe k = ", k
		!call printmat(A)
		if (A(k, k) <= ZERO) then
			!print *, "not pos def"
			! A is not positive definite!
			status = k

			! delta = np.sum(L[k, :k] ** 2) - A[k, k]
			! v = np.zeros((n,))
			! v[k] = 1.0
			! for j in range(k-1, -1, -1):  # j = k-1, k-2, ..., 0
			!     v[j] = -np.sum(L[j + 1:k+1, j] * v[j + 1:k+1]) / L[j, j]
			do j = n, 1, -1
				if (j > k) then
					v(j) = ZERO
				else if (j == k) then
					v(j) = ONE
				else
					v(j) = ZERO
					do l = j + 1, k
						v(j) = v(j) - A(l, j) * v(l)
					end do
					v(j) = v(j) / A(j, j)
				end if
			end do
			delta = -adiag(k)
			do j = 1, k - 1
				delta = delta + A(k, j) ** 2
			end do

			exit
		else
			!print *, "cholsafe, k =", k
			!call printmat(A)
			Lkk = A(k, k)
			sqrt_Lkk = sqrt(Lkk)
			do j = k + 1, n
				! A[j:, j] = A[j:, j] - A[j:, k] * A[j, k] / A[k, k]
				A(j:n, j) = A(j:n, j) - A(j:n, k) * A(j, k) / Lkk
				!do l = j, n 
				!	A(l, j) = A(l, j) - A(l, k) * A(j, k) / Lkk
				!end do
			end do
			! A[k:, k] = A[k:, k] / sqrt(A[k, k]);
			! Need to define sqrt_Lkk early since L(k,k) is updated in this loop
			A(k:n, k) = A(k:n, k) / sqrt_Lkk
			!do l = k, n 
			!	A(l, k) = A(l, k) / sqrt_Lkk
			!end do
		end if
	end do
end subroutine cholsafe

subroutine cholsolve(A, x)
	! --------------------------------------------------------------- !
	! Solve using Cholesky factorization: x <-- A \ x
	! 
	! Since A is symmetric, no difference between A \ x and A^T \ x
	! --------------------------------------------------------------- !
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : issymmetric

	implicit none

	! Inputs/outputs
	real(RP), intent(in) :: A(:,:)  ! A(N,N)
	real(RP), intent(inout) :: x(:)  ! X(N)

	! Locals
	character(len=*), parameter :: srname = 'CHOLSOLVE'
	integer(IK) :: j, k
	integer(IK) :: n

	! Sizes.
    n = int(size(A, 1), kind(n))

	! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(A, 1) == n .and. issymmetric(A), 'A is n-by-n and symmetric', srname)
		call assert(size(x) == n, 'SIZE(X) == N', srname)
    end if

	! solve L * y = b
	x(1) = x(1) / A(1,1)
	do j = 2, n 
		do k = j-1, 1, -1
			x(j) = x(j) - x(k) * A(j,k)
		end do
		x(j) = x(j) / A(j,j)
	end do
	
	! solve L^T * x = y
	x(n) = x(n) / A(n, n)
	do j=n-1, 1, -1
		do k=j+1, n 
			x(j) = x(j) - x(k) * A(k,j)
		end do		
		x(j) = x(j) / A(j,j)
	end do
	
end subroutine cholsolve

subroutine pi3(H_plus_lambda_I, x_lambda, pi, d1pi, d2pi, d3pi)
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, HALF, TWO, ONE, DEBUGGING
	use, non_intrinsic :: debug_mod, only : assert
	use, non_intrinsic :: linalg_mod, only : inprod

    implicit none

	! Inputs
	real(RP), intent(in) :: H_plus_lambda_I(:, :)  ! H_p_l_I(N,N)
	real(RP), intent(in) :: x_lambda(:)  ! x_lambda(N)
	
	! Outputs
	real(RP), intent(out) :: pi
	real(RP), intent(out) :: d1pi
	real(RP), intent(out) :: d2pi
	real(RP), intent(out) :: d3pi

	! Locals
	character(len=*), parameter :: srname = 'PI3'
	integer(IK) :: n
	real(RP), parameter :: alpha0 = ONE
	real(RP), parameter :: alpha1 = 6.0
    real(RP) :: x1(size(x_lambda))
	real(RP) :: x2(size(x_lambda))

    ! Sizes.
    n = int(size(x_lambda), kind(n))

	! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(H_plus_lambda_I, 1) == n .and. size(H_plus_lambda_I, 2) == n, 'H_p_l_I is n-by-n', srname)
		call assert(size(x_lambda) == n, 'SIZE(X) == N', srname)
    end if

	! zero-th deriv is just pi(lambda) = ||x(lambda)||^2
	pi = sum(x_lambda**2) 

	! x1 <-- (H+lambda*I) \ x1 = -1 * (H+lambda*I) \ x_lambda
	x1 = -x_lambda
	call cholsolve(H_plus_lambda_I, x1)

	! x2 <-- (H+lambda*I) \ x2 = -2 * (H+lambda*I) \ x1
	x2 = -TWO * x1
	call cholsolve(H_plus_lambda_I, x2)

	! pi'(lambda) = 2*alpha0 * dot(x_lambda, x1)
	d1pi = TWO * alpha0 * inprod(x_lambda, x1)
	! pi''(lambda) = alpha1 * ||x1||^2
	d2pi = alpha1 * sum(x1**2) 
	! pi'''(lambda) = 2*alpha1 * dot(x1, x2)
	d3pi = TWO * alpha1 * inprod(x1, x2) 
endsubroutine pi3

subroutine pi3beta(H_plus_lambda_I, x_lambda, beta, pi_beta, d1pi_beta, d2pi_beta, d3pi_beta)
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, HALF, TWO, ONE

    implicit none

	! Inputs
	real(RP), intent(in) :: H_plus_lambda_I(:, :)  ! H_p_l_I(N,N)
	real(RP), intent(in) :: x_lambda(:)  ! x_lambda(N)
	real(RP), intent(in) :: beta
	
	! Outputs
	real(RP), intent(out) :: pi_beta
	real(RP), intent(out) :: d1pi_beta
	real(RP), intent(out) :: d2pi_beta
	real(RP), intent(out) :: d3pi_beta

	! Locals
	real(RP) :: pi, d1pi, d2pi, d3pi
	real(RP) :: half_beta

	half_beta = HALF * beta

	call pi3(H_plus_lambda_I, x_lambda, pi, d1pi, d2pi, d3pi)

	pi_beta = pi ** half_beta
	d1pi_beta = half_beta * (pi ** (half_beta - ONE)) * d1pi
	d2pi_beta = half_beta * (pi ** (half_beta - ONE)) * d2pi + half_beta * (half_beta - ONE) * (pi**(half_beta-TWO)) * (d1pi**2)
	d3pi_beta = (pi**2) * d3pi + 3.0 * (half_beta - ONE) * pi * d1pi * d2pi + (half_beta - ONE) * (half_beta - TWO) * (d1pi**3)
	d3pi_beta = d3pi_beta * half_beta * (pi**(half_beta-3.0))

end subroutine pi3beta

subroutine newlda(H_plus_lambda_I, x_lambda, lambdaC, delta, beta, k, new_lambda, is_tr, status)
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HALF, ONE, TWO, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert

	implicit none

	! Inputs
	real(RP), intent(in) :: H_plus_lambda_I(:, :)  ! H_p_l_I(N,N)
	real(RP), intent(in) :: x_lambda(:)  ! X_LAMBDA(N)
	real(RP), intent(in) :: lambdaC
	real(RP), intent(in) :: delta
	real(RP), intent(in) :: beta
	integer(IK), intent(in) :: k
	logical, intent(in) :: is_tr

	! Outputs
	real(RP), intent(out) :: new_lambda
	integer(IK), intent(out) :: status

	! Locals
	real(RP) :: pi_beta, d1pi_beta, d2pi_beta, d3pi_beta
	real(RP) :: sq_delta
	real(RP) :: d1, d2, d3
	real(RP) :: p0, p1, p2, p3
	integer(IK) :: nroots

	! Check for valid combinations of (beta,k)
	status = 0
	if (k < 1 .or. k > 3 .or. beta == ZERO) then
		status = -1
	else if (.not. is_tr) then
		! ARC only implements specific combinations
		if (.not. (beta == -ONE .and. k == 1) .and. .not. (beta == TWO .and. k == 2) .and. .not. (beta == TWO .and. k == 3)) then
			status = -1
		end if
	end if
	
	if (status == 0) then
		sq_delta = delta * delta
		call pi3beta(H_plus_lambda_I, x_lambda, beta, pi_beta, d1pi_beta, d2pi_beta, d3pi_beta)

		if (k == 1) then
			! For (TRS), solve: pi_beta + d1pi_beta*d - delta^beta = 0
			!
			! For (ARC), then have beta=-1, and solve
			!    pi_beta + d1pi_beta*d - (lambdaC+d)^(-1) / delta^(-1) = 0
			! or
			!    d1pi_beta * d^2 + (pi_beta + d1pi_beta*lambdaC) * d + (lambdaC*pi_beta - delta) = 0
			if (is_tr) then
				call cubicroots(ZERO, ZERO, d1pi_beta, pi_beta - (delta ** beta), d1, d2, d3, nroots)
			else
				call cubicroots(ZERO, d1pi_beta, pi_beta + d1pi_beta * lambdaC, lambdaC * pi_beta - delta, d1, d2, d3, nroots)
			end if
		else if (k == 2) then
			! For (TRS), solve: pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 - delta^beta = 0
			!
			! For (ARC), then have beta=2, so solve:
			!   pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 - (lambdaC+d)^2 / delta^2 = 0
			! or
			!   (pi_beta - lambdaC^2 / delta^2) + (d1pi_beta - 2*lambdaC/delta^2) * d + (0.5*d2pi_beta - 1/delta^2) * d^2 = 0
			if (is_tr) then
				call cubicroots(ZERO, HALF * d2pi_beta, d1pi_beta, pi_beta - (delta ** beta), d1, d2, d3, nroots)
			else
				p0 = ZERO
				p1 = HALF * d2pi_beta - ONE / (sq_delta)
				p2 = d1pi_beta - TWO * lambdaC / sq_delta
				p3 = pi_beta - lambdaC * lambdaC / sq_delta
				call cubicroots(p0, p1, p2, p3, d1, d2, d3, nroots)
			end if
		else
			! For (TRS), solve: pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 + (1/6) * d3pi_beta * d^3 - delta^beta = 0
			!
			! For (ARC), then have beta=2, so solve:
			!     pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 + (1/6) * d3pi_beta * d^3 - (lambdaC+d)^2 / delta^2 = 0
			! or
			!     (pi_beta - lambdaC^2 / delta^2) + (d1pi_beta - 2*lambdaC/delta^2) * d + (0.5*d2pi_beta - 1/delta^2) * d^2 + (1/6) * d3pi_beta * d^3 = 0
			if (is_tr) then
				call cubicroots(d3pi_beta / 6.0, HALF * d2pi_beta, d1pi_beta, pi_beta - (delta ** beta), d1, d2, d3, nroots)
			else
				p0 = d3pi_beta / 6.0
				p1 = HALF * d2pi_beta - ONE / (sq_delta)
				p2 = d1pi_beta - TWO * lambdaC / sq_delta
				p3 = pi_beta - lambdaC * lambdaC / sq_delta
				call cubicroots(p0, p1, p2, p3, d1, d2, d3, nroots)
			end if
		end if

		if (nroots == 0) then
			new_lambda = lambdaC
			status = -1
		else 
			status = 0
			if (nroots == 1) then
				new_lambda = lambdaC + d1
			else if (nroots == 2) then
				new_lambda = lambdaC + d2
			else
				new_lambda = lambdaC + d3
			end if
		end if
	end if
	
end subroutine newlda

subroutine hardstep(xs, ztmp, delta, lambdaC, alpha, is_tr, status)
	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, TWO, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
	use, non_intrinsic :: linalg_mod, only : inprod

	implicit none

	! Inputs
	real(RP), intent(in) :: xs(:)  ! XS(N)
	real(RP), intent(in) :: ztmp(:) ! ZTMP(N)
	real(RP), intent(in) :: delta
	real(RP), intent(in) :: lambdaC
	logical, intent(in) :: is_tr

	! Outputs
	real(RP), intent(out) :: alpha
	integer(IK), intent(out) :: status

	! Locals
	character(len=*), parameter :: srname = 'HARDSTEP'
	integer(IK) :: n
	integer(IK) :: nroots
	real(RP) :: alpha1, alpha2, rhs

	! Sizes.
    n = int(size(xs), kind(n))

    ! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(xs) == n, 'SIZE(XS) == N', srname)
        call assert(size(ztmp) == n, 'SIZE(ZTMP) == N', srname)
    end if
	
	if (is_tr) then
		rhs = delta
	else
		rhs = lambdaC / delta
	end if
	
	call quadroots(sum(ztmp**2), TWO * inprod(xs, ztmp), sum(xs**2) - rhs * rhs, alpha1, alpha2, nroots)

	if (nroots > 0) then
		! any root is fine, since sign of z is arbitrary
		alpha = alpha1  
		status = 0
	else
		status = -1
	end if

end subroutine hardstep

subroutine cauchy(hess, g, delta, is_tr, s)
    ! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HALF, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : matprod, inprod, issymmetric

    implicit none

    ! Inputs
    real(RP), intent(in):: hess(:, :)  ! HESS(N,N)
    real(RP), intent(in):: g(:)  ! G(N,N)
    real(RP), intent(in):: delta
    logical, intent(in) :: is_tr

    ! Outputs
    real(RP), intent(out) :: s(:)  ! S(N)

    ! Local variables
    character(len=*), parameter :: srname = 'CAUCHY'
    integer(IK) :: n
    real(RP) :: alpha
    real(RP) :: gHg
    real(RP) :: normg
    real(RP) :: tmp
    real(RP) :: hg(size(g))
	real(RP) :: tau

    ! Sizes.
    n = int(size(g), kind(n))

    ! Preconditions
    if (DEBUGGING) then
        call assert(n >= 1, 'N >= 1', srname)
        call assert(size(g) == n, 'SIZE(G) == N', srname)
        call assert(size(hess, 1) == n .and. issymmetric(hess), 'HESS is n-by-n and symmetric', srname)
        call assert(size(s) == n, 'SIZE(S) == N', srname)
    end if

    normg = sqrt(sum(g**2))

    if (normg == ZERO) then
        ! arbitrary, since multiplying by zero vector
        alpha = ZERO
    else
        hg = matprod(hess, g)
        gHg = inprod(g, hg)
        tmp = gHg / normg**3
        if (is_tr) then
            if (gHg <= ZERO) then
                tau = ONE
            else
                tau = min(ONE / (delta * tmp), ONE)
            end if
            alpha = -tau * delta / normg
        else
            alpha = -HALF * tmp / delta + HALF * sqrt(tmp * tmp + 4 * delta / normg) / delta
        end if
    end if

    s = alpha * g
end subroutine cauchy

subroutine quadroots(p0, p1, p2, x1, x2, nroots)
    ! ------------------------------------------------------------------------------- !
    ! Find the real roots x1 <= x2 to the quadratic equation p0*x^2 + p1*x + p2 = 0,
	! and return the number of real roots.
	! 
	! If 1 real root, it is given in both x1 and x2
	! 
	! If all p0=p1=p2=0, then returns 0 and leaves x1,x2 unchanged.
	! 
	! Based on function QUADROOTS, from
	! T. R. F. Nonweiler. Algorithm 326: Roots of Low-Order Polynomial Equations
	! Communications of the ACM, 11:4 (1968), pp. 269-270.
	! https://dl.acm.org/doi/10.1145/362991.363039
	! 
	! Computed solutions are refined using 1 iteration of Newton's method 
    ! (an idea used in GALAHAD/roots.f90)
    ! ------------------------------------------------------------------------------- !

    ! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, HALF, ZERO

    implicit none

    ! Inputs
    real(RP), intent(in) :: p0
    real(RP), intent(in) :: p1
    real(RP), intent(in) :: p2

    ! Outputs
    real(RP), intent(out) :: x1
    real(RP), intent(out) :: x2
    integer(IK), intent(out) :: nroots

    ! Localc
    real(RP) :: b, c, d
    real(RP) :: p, dp

	! Some manual overrides
	if (p0 == ZERO) then
		if (p1 == ZERO) then
			! solve p2=0, either zero or infinitely many solutions, but return no solutions in this case
			nroots = 0
		else
			! Linear equation p1*x + p2 = 0
			x1 = -p2 / p1
			x2 = x1
			nroots = 1
        end if
	else 
        if (p2 == ZERO) then
		    ! p0*x^2 + p1*x = 0
            if (p1 == ZERO) then
                x1 = ZERO;
                x2 = x1
                nroots = 1
            else
                x1 = ZERO
                x2 = -p1 / p0
                nroots = 2
            end if
        else  
            ! p0 and p2 are both nonzero, the usual case
            b = -HALF * p1 / p0
            c = p2 / p0
            d = b * b - c
            if (d > ZERO) then
                ! 2 real roots
                if (b > ZERO) then
                    x1 = b + sqrt(d)
                else
                    x1 = b - sqrt(d)
				end if
                x2 = c / x1
                nroots = 2
            else 
                if (d == ZERO) then
                    ! 1 real root
                    x1 = b
                    x2 = x1
                    nroots = 1
                else
                    ! no real roots, don't set x1 or x2
                    nroots = 0
                end if
            end if
        end if
    end if

	! Put in increasing order, x1 <= x2
	if (nroots >= 2 .and. x1 > x2) then
        ! Swap x1 and x2 (using b as tmp variable)
        b = x2
		x2 = x1
		x1 = b
    end if

	! 1 iteration of Newton's method to refine each root
	! p(x) = p0*x^2 + p1*x + p2
	! p'(x) = 2*p0*x + p1
	if (nroots >= 1) then
		! refine x1
		p = p0 * x1 * x1 + p1 * x1 + p2;
		dp = 2 * p0 * x1 + p1;
		if (dp .ne. ZERO) then
            x1 = x1 - p / dp;
        end if
    end if

	if (nroots >= 2) then
		! refine x2
		p = p0 * x2 * x2 + p1 * x2 + p2;
		dp = 2 * p0 * x2 + p1;
		if (dp .ne. ZERO) then
            x2 = x2 - p / dp;
        end if
    end if
end subroutine quadroots


subroutine cubicroots(p0, p1, p2, p3, x1, x2, x3, nroots)
	! ------------------------------------------------------------------------------------------- !
	! Find the real roots x1 <= x2 <= x3 to the cubic equation p0*x^3 + p1*x^2 + p2*x + p3 = 0,
	! and return the number of real roots.
	! 
	! If any repeated roots are found, these are all returned.
	! 
	! If all p0=p1=p2=p3=0, then returns 0 and leaves x1,x2,x3 unchanged.
	! 
	! Based on function CUBICROOTS, from
	! T. R. F. Nonweiler. Algorithm 326: Roots of Low-Order Polynomial Equations
	! Communications of the ACM, 11:4 (1968), pp. 269-270.
	! https://dl.acm.org/doi/10.1145/362991.363039
	! 
	! Computed solutions are refined using 1 iteration of Newton's method 
	! (an idea used in GALAHAD/roots.f90)
	! ------------------------------------------------------------------------------------------- !

    ! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, HALF, ONE, TWO, ZERO

    implicit none

    ! Inputs
    real(RP), intent(in) :: p0
    real(RP), intent(in) :: p1
    real(RP), intent(in) :: p2
    real(RP), intent(in) :: p3

    ! Outputs
    real(RP), intent(out) :: x1
    real(RP), intent(out) :: x2
    real(RP), intent(out) :: x3
    integer(IK), intent(out) :: nroots

    ! Locals
    real(RP) :: q1, q2, q3
    real(RP) :: s, t, b, c, d
    real(RP) :: p, dp

	if (p0 == ZERO) then
		call quadroots(p1, p2, p3, x1, x2, nroots)
	else if (p3 == ZERO) then
		! p(x) = p0*x^3 + p1*x^2 + p2*x --> x=0 is one root
		x1 = 0;
		call quadroots(p0, p1, p2, x2, x3, nroots)
        nroots = nroots + 1
	else
		! Get normalized equation: x^3 + q1*x^2 + q2*x + q3 = 0
		q1 = p1 / p0
		q2 = p2 / p0
		q3 = p3 / p0

		s = q1 / 3.0
		t = s * q1
		b = HALF * (s * (t / 1.5 - q2) + q3)
		t = (t - q2) / 3.0
		c = t * t * t
		d = b * b - c

		if (d >= ZERO) then
			! 1 real root + other 2 roots are either double real root or complex conjugate pair 
			d = (sqrt(d) + abs(b)) ** (ONE / 3.0)
			if (d .ne. ZERO) then
                if (b > ZERO) then
                    b = -d 
                else
                    b = d
                end if
				c = t / b
            end if
			d = sqrt(0.75) * (b - c)
			b = b + c
			c = -HALF * b - s
		
			! roots are x=b-s and x=c +/- d*i (real double root if d=0)
			x1 = b - s
			if (d == ZERO) then
				x2 = c
				x3 = x2
				nroots = 3
			else
				nroots = 1
            end if
		else
			! 3 distinct roots
            if (b == ZERO) then
                d = atan(ONE) / 1.5
            else
                d = atan(sqrt(-d) / abs(b)) / 3.0
            end if
            if (b < ZERO) then
                b = TWO * sqrt(t)
            else
                b = -TWO * sqrt(t)
            end if
			c = cos(d) * b
			t = -sqrt(0.75) * sin(d) * b - HALF * c
			d = -t - c - s
			c = c - s
			t = t - s
			x1 = c
			x2 = d
			x3 = t
			nroots = 3
        end if
    end if

	! Make sure roots are sorted in increasing order, x1 <= x2 <= x3 (where relevant)
	if (nroots > 1) then
		if (x1 > x2) then
			b = x2
			x2 = x1
			x1 = b
        end if
		if ((nroots > 2) .and. (x2 > x3)) then
			b = x3
			x3 = x2
			x2 = b

			! if x3 (now x2) was originally the smallest, now might need to do another swap with x1
			if (x1 > x2) then
				b = x2
				x2 = x1
				x1 = b
            end if
        end if
    end if

	! 1 iteration of Newton's method to refine each root
	! p(x) = p0*x^3 + p1*x^2 + p2*x + p3
	! p'(x) = 3*p0*x^2 + 2*p1*x + p2
	
	if (nroots >= 1) then
		! refine x1
		p = p0 * x1 * x1 * x1 + p1 * x1 * x1 + p2 * x1 + p3
		dp = 3 * p0 * x1 * x1 + 2 * p1 * x1 + p2
		if (dp .ne. ZERO) then
            x1 = x1 - p / dp
        end if
    end if

	if (nroots >= 2) then
		! refine x2
		p = p0 * x2 * x2 * x2 + p1 * x2 * x2 + p2 * x2 + p3
		dp = 3 * p0 * x2 * x2 + 2 * p1 * x2 + p2
		if (dp .ne. ZERO) then
            x2 = x2 - p / dp
        end if
    end if

	if (nroots >= 3) then
		! refine x3
		p = p0 * x3 * x3 * x3 + p1 * x3 * x3 + p2 * x3 + p3
		dp = 3 * p0 * x3 * x3 + 2 * p1 * x3 + p2
		if (dp .ne. ZERO) then
            x3 = x3 - p / dp
        end if
    end if

end subroutine cubicroots

subroutine testquad()
	! Tests for quadroots

	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert

	implicit none

	real(RP), parameter :: thresh = 1e-10
	real(RP) :: p0, p1, p2, x1, x2
	real(RP) :: x1true, x2true
	integer(IK) :: nroots, nroots_true, nroots_true_with_repeats
	real(RP) :: xerr, ferr

	print *, "************** TESTQUAD **************"

	print *, "Basic1"
	p0 = 2
	p1 = 0
	p2 = -1
	x1true = -0.707106781186548
	x2true = 0.707106781186548
	nroots_true = 2
	nroots_true_with_repeats = 2
	call quadroots(p0, p1, p2, x1, x2, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**2 + p1 * x1 + p2)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**2 + p1 * x2 + p2)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "Double Root"
	p0 = 1
	p1 = 0
	p2 = 0
	x1true = 0
	x2true = 0
	nroots_true = 1
	nroots_true_with_repeats = 2
	call quadroots(p0, p1, p2, x1, x2, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**2 + p1 * x1 + p2)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**2 + p1 * x2 + p2)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "Linear"
	p0 = 0
	p1 = 2
	p2 = 1
	x1true = -0.5
	x2true = 1E10
	nroots_true = 1
	nroots_true_with_repeats = 1
	call quadroots(p0, p1, p2, x1, x2, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**2 + p1 * x1 + p2)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**2 + p1 * x2 + p2)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "Constant"
	p0 = 0
	p1 = 0
	p2 = 1
	x1true = 1E10
	x2true = 1E10
	nroots_true = 0
	nroots_true_with_repeats = 0
	call quadroots(p0, p1, p2, x1, x2, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**2 + p1 * x1 + p2)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**2 + p1 * x2 + p2)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "Basic2"
	p0 = 2
	p1 = 0
	p2 = -2
	x1true = -1
	x2true = 1
	nroots_true = 2
	nroots_true_with_repeats = 2
	call quadroots(p0, p1, p2, x1, x2, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**2 + p1 * x1 + p2)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**2 + p1 * x2 + p2)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "Basic3"
	p0 = 3
	p1 = 6
	p2 = -9
	x1true = -3
	x2true = 1
	nroots_true = 2
	nroots_true_with_repeats = 2
	call quadroots(p0, p1, p2, x1, x2, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**2 + p1 * x1 + p2)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**2 + p1 * x2 + p2)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "************** END TESTQUAD **************"
end subroutine testquad

subroutine testcubic()
	! Tests for cubicroots

	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert

	implicit none

	real(RP) :: p0, p1, p2, p3, x1, x2, x3
	real(RP) :: x1true, x2true, x3true
	integer(IK) :: nroots, nroots_true, nroots_true_with_repeats
	real(RP) :: xerr, ferr

	print *, "************** TESTCUBIC **************"

	print *, "Basic1"
	p0 = 1
	p1 = -6
	p2 = 4
	p3 = 12
	x1true = -1.05137424173104
	x2true = 2.51730404500831
	x3true = 4.53407019672273
	nroots_true = 3
	nroots_true_with_repeats = 3
	call cubicroots(p0, p1, p2, p3, x1, x2, x3, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**3 + p1 * x1 ** 2 + p2 * x1 + p3)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**3 + p1 * x2**2 + p2 * x2 + p3)
	end if
	if (nroots >= 3) then
		xerr = max(xerr, x3 - x3true)
		ferr = max(ferr, p0 * x3**3 + p1 * x3 ** 2 + p2 * x3 + p3)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "Basic2"
	p0 = 3
	p1 = -7.5
	p2 = -16.5
	p3 = 21
	x1true = -2
	x2true = 1
	x3true = 3.5
	nroots_true = 3
	nroots_true_with_repeats = 3
	call cubicroots(p0, p1, p2, p3, x1, x2, x3, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**3 + p1 * x1 ** 2 + p2 * x1 + p3)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**3 + p1 * x2**2 + p2 * x2 + p3)
	end if
	if (nroots >= 3) then
		xerr = max(xerr, x3 - x3true)
		ferr = max(ferr, p0 * x3**3 + p1 * x3 ** 2 + p2 * x3 + p3)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "TripleRoot"
	p0 = 1
	p1 = 0
	p2 = 0
	p3 = 0
	x1true = 0
	x2true = 0
	x3true = 0
	nroots_true = 1
	nroots_true_with_repeats = 3
	call cubicroots(p0, p1, p2, p3, x1, x2, x3, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**3 + p1 * x1 ** 2 + p2 * x1 + p3)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**3 + p1 * x2**2 + p2 * x2 + p3)
	end if
	if (nroots >= 3) then
		xerr = max(xerr, x3 - x3true)
		ferr = max(ferr, p0 * x3**3 + p1 * x3 ** 2 + p2 * x3 + p3)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "SingleRoot"
	p0 = 1
	p1 = 0
	p2 = 0
	p3 = 2
	x1true = -1.25992104989487
	x2true = 0
	x3true = 0
	nroots_true = 1
	nroots_true_with_repeats = 1
	call cubicroots(p0, p1, p2, p3, x1, x2, x3, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**3 + p1 * x1 ** 2 + p2 * x1 + p3)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**3 + p1 * x2**2 + p2 * x2 + p3)
	end if
	if (nroots >= 3) then
		xerr = max(xerr, x3 - x3true)
		ferr = max(ferr, p0 * x3**3 + p1 * x3 ** 2 + p2 * x3 + p3)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "DoubleRoot"
	p0 = 1
	p1 = -1
	p2 = -1
	p3 = 1
	x1true = -1
	x2true = 1
	x3true = 1
	nroots_true = 2
	nroots_true_with_repeats = 3
	call cubicroots(p0, p1, p2, p3, x1, x2, x3, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**3 + p1 * x1 ** 2 + p2 * x1 + p3)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**3 + p1 * x2**2 + p2 * x2 + p3)
	end if
	if (nroots >= 3) then
		xerr = max(xerr, x3 - x3true)
		ferr = max(ferr, p0 * x3**3 + p1 * x3 ** 2 + p2 * x3 + p3)
	end if
	print *, "xerr =", xerr, "ferr =", ferr

	print *, "Quadratic"
	p0 = 0
	p1 = 2
	p2 = 0
	p3 = -1
	x1true = -0.707106781186548
	x2true = 0.707106781186548
	x3true = 0
	nroots_true = 2
	nroots_true_with_repeats = 2
	call cubicroots(p0, p1, p2, p3, x1, x2, x3, nroots)
	if (nroots < nroots_true .or. nroots > nroots_true_with_repeats) then
		print *, "bad nroots = ", nroots
	end if
	xerr = ZERO
	ferr = ZERO
	if (nroots >= 1) then
		xerr = max(xerr, x1 - x1true)
		ferr = max(ferr, p0 * x1**3 + p1 * x1 ** 2 + p2 * x1 + p3)
	end if
	if (nroots >= 2) then
		xerr = max(xerr, x2 - x2true)
		ferr = max(ferr, p0 * x2**3 + p1 * x2**2 + p2 * x2 + p3)
	end if
	if (nroots >= 3) then
		xerr = max(xerr, x3 - x3true)
		ferr = max(ferr, p0 * x3**3 + p1 * x3 ** 2 + p2 * x3 + p3)
	end if
	print *, "xerr =", xerr, "ferr =", ferr


	print *, "************** END TESTCUBIC **************"
end subroutine testcubic

subroutine testchol()
	! Tests for Cholesky routines (cholsafe, cholsolve)

	! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : issymmetric

	implicit none

	integer(IK), parameter :: n = 3
	real(RP), parameter :: thresh = 1.0E-7

	real(RP) :: A(n, n), L(n, n), Aorig(n, n), b(n), borig(n), v(n), vtrue(n), xtrue(n)
	real(RP) :: delta, delta_true, err
	integer(IK) :: status
	integer(IK) :: i, j

	print *, "************** TESTCHOL **************"

	A = reshape([10.5, 2.5, 3.5, 2.5, 50.0, 6.0, 3.5, 6.0, 90.0], [n,n])
	Aorig = A
	L = transpose(reshape([3.24037035, 2.5, 3.5, 0.77151675, 7.0288521, 6.0, 1.08012345, 0.7350655, 9.39643614], [n,n]))

	delta_true = ZERO
	vtrue = ZERO
	call cholsafe(A, v, delta, status)
	print *, "CHOLSAFE - case 1"
	print *, "status (expect 0) = ", status
	!print *, "L = "
	!call printmat(A)
	!print *, "Ltrue = "
	!call printmat(L)
	print *, "L err = ", maxval(abs(A - L))
	print *, "delta err = ", abs(delta - delta_true)
	print *, "v err = ", maxval(abs(v - vtrue))

	print *, "CHOLSOLVE - case 1"
	b = (/ 1.0, -1.0, 3.0 /)
	xtrue = (/ 8376.0 / 91604.0, -2599.0 / 91604.0, 2901.0 / 91604.0 /)
	call cholsolve(A, b)
	print *, "x err = ", maxval(abs(b - xtrue))

	print *, "CHOLSAFE - case 2"
	A = Aorig
	A(2, 2) = -50.0
	vtrue = (/ -0.23809524, 1.0, 0.0 /)
	delta_true = 50.595238095238095
	!call printmat(A)
	call cholsafe(A, v, delta, status)
	print *, "status (expect 2) = ", status
	print *, "delta err = ", abs(delta - delta_true)
	print *, "v err = ", maxval(abs(v - vtrue))

	print *, "************** END TESTCHOL **************"
end subroutine testchol

subroutine trsunc(delta, g_in, hess_in, lambda, s)
    use, non_intrinsic :: consts_mod, only : RP

    implicit none

    ! Inputs
    real(RP), intent(in) :: delta
    real(RP), intent(in) :: g_in(:)  ! G_IN(N)
    real(RP), intent(in) :: hess_in(:, :)  ! HESS_IN(N, N)

    ! Outputs
    real(RP), intent(out) :: lambda
    real(RP), intent(out) :: s(:)  ! S(N)

    call trglob(delta, g_in, hess_in, lambda, s, .true.)

	!call testchol()
	!call testquad()
	!call testcubic()
end subroutine trsunc

subroutine arcunc(delta, g_in, hess_in, lambda, s)
    use, non_intrinsic :: consts_mod, only : RP

    implicit none

    ! Inputs
    real(RP), intent(in) :: delta
    real(RP), intent(in) :: g_in(:)  ! G_IN(N)
    real(RP), intent(in) :: hess_in(:, :)  ! HESS_IN(N, N)

    ! Outputs
    real(RP), intent(out) :: lambda
    real(RP), intent(out) :: s(:)  ! S(N)

    call trglob(delta, g_in, hess_in, lambda, s, .false.)
end subroutine arcunc

end submodule trsunc_mod