submodule (trustregion_mod) trsunc_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Thursday, April 04, 2024 PM09:26:23
!--------------------------------------------------------------------------------------------------!

implicit none

contains

subroutine trglob(delta, g_in, hess_in, lambda, s, is_tr)
! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, HALF, REALMIN, ZERO, TENTH, EPS, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : inprod, issymmetric, norm, project, matprod
use, non_intrinsic :: univar_mod, only : circle_min

use, non_intrinsic :: qdec_mod, only : qdec, crdec

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: g_in(:)  ! G_IN(N)
real(RP), intent(in) :: hess_in(:, :)  ! HESS_IN(N, N)
logical, intent(in) :: is_tr

! Outputs
real(RP), intent(out) :: lambda
real(RP), intent(out) :: s(:)  ! S(N)

! Locals
character(len=*), parameter :: srname = 'TRGLOB'
integer(IP) :: n
logical :: scaled
real(RP) :: gopt(size(g_in))
real(RP) :: hess(size(hess_in, 1), size(hess_in, 2))
real(RP) :: modscal

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

if (is_tr .and. delta < boundary_thresh) then 
    ! Easy quit: if constraint is ||s|| <= 0, then only feasible solution is s=0
    s = 0.0
    lambda = 0.0
    return
end if

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

	bool result = false;  ! have we found a value for x yet?

	! Find initial lambda region and current guess
	Number lambdaL, lambdaU, lambdaC;
	Index status = 0;

	initial_lambda_region(H_scal, g_scal, delta, lambdaL, lambdaU, is_tr);
	lambdaC = (lambdaL == 0.0 ? 0.0 : std::max(gamma * sqrt(lambdaL * lambdaU), lambdaL + theta * (lambdaU - lambdaL)));
	DFOPT_LOG_VERBOSE_IF(debug) << "Initially, lambdaL = " << lambdaL << " and lambdaU = " << lambdaU;

	! Main loop
	Number dval; ! output from safe_cholesky
	Number norm_check;  ! desired value of ||x(lambda)||, potentially updated at each iteration

	! For (TRS), check for interior solution before proceeding further
	if (is_tr && lambdaL == 0.0)
	{
		status = solve_kkt_system(H_scal, g_scal, 0.0, x, dval);
		if (status == 0 && x.norm2() <= delta + boundary_thresh * std::max(1.0, delta))
		{
			DFOPT_LOG_VERBOSE_IF(debug) << "Interior solution found, terminating";
			return true;
		}
	}

	! Initialize z to a unit vector orthogonal to g (used for all inverse iterations)
	Number rho;
	Number inv_sqrt_n = 1.0 / sqrt((Number)n);
	z.fill_with(inv_sqrt_n);  ! z[:] = inv_sqrt_n

	! First, make sure z is not parallel to g...
	Number gnorm = g_scal.norm2();
	if (gnorm != 0.0)
	{
		! If g=0, then z is definitely not parallel to g
		Number cos_z_g = dot(z, g_scal) / g_scal.norm2();
		if (std::abs(cos_z_g) >= 1 - z_parallel_g_thresh)
			! cos(z,g) ~ 1, so g is basically a vector with all entries the same...
			z[0] = -inv_sqrt_n;  ! this will definitely make z not parallel to g
		! Now make z orthogonal to g using the formula z = z - (ghat.T * z) * ghat, where ghat = g / ||g||
		z.add_multiple(-dot(z, g_scal) / g_scal.sqnorm2(), g_scal);  ! z += [-dot(g,z)/||g||^2] * g
	}



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

	bool found_lambdaC_in_L = false;
	bool potential_hard_case = false;
	bool hard_case = false;
	Index it = 0;

	while (std::abs(lambdaU - lambdaL) > std::max(bisection_thresh * std::max(std::abs(lambdaL), std::abs(lambdaU)), bisection_thresh)
		&& (!found_lambdaC_in_L) && (!hard_case))
	{
		if (lambdaL >= lambdaU) break;  ! always expect lambdaL < lambdaU, so stop the loop if this is violated

		DFOPT_LOG_VERBOSE_IF(debug) << "It " << it << ", (lambdaL, lambdaC, lambdaU) = " << lambdaL << ", " << lambdaC << ", " << lambdaU;
		it++;

		if (!potential_hard_case)
		{
			! Identify which region lambdaC is in
			status = solve_kkt_system(H_scal, g_scal, lambdaC, x, dval);

			if (status == 0)
			{
				! H + lambdaC*I is positive definite
				if (x.norm2() >= (is_tr ? delta : lambdaC / delta))
				{
					! lambdaC in L, i.e. ||x(lambda)|| >= delta or lambdaC/delta
					! This is the good region, where we switch to the fast-converging phase 2
					DFOPT_LOG_VERBOSE_IF(debug) << "lambdaC in L, stopping";
					found_lambdaC_in_L = true;
					break; ! exit loop
				}
				else
				{
					! lambda in G, i.e. ||x(lambda)|| < delta (lambda too large)
					lambdaU = std::min(lambdaU, lambdaC);

					! Update lambdaL using inverse iteration
					for (Index i = 0; i < num_inverse_iters_lambda_in_G; ++i)
					{
						! H_plus_lambda_I has already set and factorized in solve_kkt_system above
						H_plus_lambda_I.cholesky_solve(z);  ! z <-- (H+lambdaC*I) \ z
						z.scale(1.0 / z.norm2());  ! z <-- z / ||z||
					}

					rho = H_scal.quadform(z);
					lambdaL = std::max(lambdaL, -rho);

					! Get under-approximations of lambdaC using Taylor approximations, and update lambdaL based on this
					Number lambda1_neg1, lambda2_2, lambda3_2, lambdaT;
					lambdaT = lambdaL;
					status = taylor_new_lambda(x, lambdaC, delta, -1.0, 1, lambda1_neg1, is_tr);
					if (status == 0) lambdaT = std::max(lambdaT, lambda1_neg1);
					status = taylor_new_lambda(x, lambdaC, delta, 2.0, 2, lambda2_2, is_tr);
					if (status == 0) lambdaT = std::max(lambdaT, lambda2_2);
					status = taylor_new_lambda(x, lambdaC, delta, 2.0, 3, lambda3_2, is_tr);
					if (status == 0) lambdaT = std::max(lambdaT, lambda3_2);

					Number small_width = small_width_thresh * std::abs(lambdaU - lambdaL);
					if ((lambdaT >= lambdaL + small_width) && (lambdaT <= lambdaU - small_width))
						lambdaC = lambdaT;
					else
						lambdaC = std::max(gamma * sqrt(lambdaL * lambdaU), lambdaL + theta * (lambdaU - lambdaL));
					!DFOPT_LOG_VERBOSE_IF(debug) << "lambda in G";
				}
			}
			else if (status > 0)
			{
				! lambda in N, i.e. H+lambdaC*I is not positive definite (lambdaC too small)
				! Here, x and d are set so that x.T * (H + lambdaC*I + d * ek * ek.T) * x = 0
				lambdaL = std::max(lambdaL, lambdaC);
				lambdaL = std::max(lambdaL, dval / x.sqnorm2() + lambdaC);
				lambdaC = std::max(gamma * sqrt(lambdaL * lambdaU), lambdaL + theta * (lambdaU - lambdaL));
				!DFOPT_LOG_VERBOSE_IF(debug) << "lambda in N";
			}
			else
			{
				DFOPT_LOG_VERBOSE_IF(debug) << "Error in solve_kkt_system, stopping";
				break;
			}

			! Check if we are in the potential hard case now
			if (std::abs(lambdaU - lambdaL) <= std::max(potential_hard_case_thresh * std::max(std::abs(lambdaL), std::abs(lambdaU)), potential_hard_case_thresh))
			{
				DFOPT_LOG_VERBOSE_IF(debug) << " - Identified potential hard case, lambdaL = " << lambdaL << ", lambdaU = " << lambdaU;
				potential_hard_case = true;
			}
			else {
				potential_hard_case = false;
			}
		} ! end regular bisection step
		else
		{
			DFOPT_LOG_VERBOSE_IF(debug) << " - Potential hard case";
			Number old_lambdaC = lambdaC;
			lambdaC = lambdaU; ! definitely an overestimate of -lambda1

			for (Index it2 = 0; it2 < num_potential_hard_case_iters; ++it2)
			{
				DFOPT_LOG_VERBOSE_IF(debug) << " - It " << it2 << ", trying lambdaC = " << lambdaC;
				status = solve_kkt_system(H_scal, g_scal, lambdaC, x, dval);

				if (status == 0)
				{
					! H + lambdaC*I is positive definite (expected)
					if (x.norm2() >= (is_tr ? delta : lambdaC / delta))
					{
						! lambdaC in L, i.e. ||x(lambda)|| >= delta or lambdaC/delta
						! This is the good region, where we switch to the fast-converging phase 2
						DFOPT_LOG_VERBOSE_IF(debug) << "lambdaC in L, stopping";
						found_lambdaC_in_L = true;
						break; ! exit loop
					}
					! Otherwise, update lambdaC using inverse iteration
					Index nk = (it > 5 ? 1 : 2);  ! don't need too many iterations
					for (Index i = 0; i < nk; ++i)
					{
						! H_plus_lambda_I has already set and factorized in solve_kkt_system above
						H_plus_lambda_I.cholesky_solve(z);  ! z <-- (H+lambdaC*I) \ z
						z.scale(1.0 / z.norm2());  ! z <-- z / ||z||
					}

					! Update lambdaC
					rho = H_scal.quadform(z);
					Number gammak = (nk == 1 ? 1.5 : 3);
					Number new_lambda = -rho + omega * std::pow(lambdaC + rho, gammak);
					if (std::abs(new_lambda - lambdaC) < lambda_convergence_thresh)
					{
						hard_case = true;
						break;
					}
					if (new_lambda < lambdaC)
					{
						lambdaC = new_lambda;
					}
					else
					{
						! If we start lambdaC sufficiently close to a good value, this should never happen
						! So, stop the 'potential hard case' iteration and go back to regular bisection phase
						DFOPT_LOG_VERBOSE_IF(debug) << " - Started nearly hard case too soon, trying again";
						potential_hard_case = false;
						lambdaC = old_lambdaC;
						potential_hard_case_thresh *= potential_hard_case_thresh_decrease;
						break;
					}
				}
				else if (status > 0)
				{
					! H + lambdaC*I is not positive definite
					! This happens when lambdaC is slightly smaller than -lambda_{max}(H), due to rounding errors
					DFOPT_LOG_VERBOSE_IF(debug) << "Potential hard case failure (lambdaC too small from rounding errors) - treat as hard case";
					hard_case = true;
					break;
				}
				else
				{
					DFOPT_LOG_VERBOSE_IF(debug) << "Error in solve_kkt_system, stopping";
					break;
				}
			}  ! end potential hard case iteration loop

			! Stop the bisection loop
			if (found_lambdaC_in_L || hard_case || status < 0) break;

			! If the potential hard case didn't find lambdaC in L, then lambdaC is still in G and is a potential upper bound
			lambdaU = std::min(lambdaC, lambdaU);
		} ! end potential hard case
	} ! end bisection loop

	! End of main lambdaC search phase
	! Three possibilities to get here:
	! 1. Found lambdaC in L, which needs to be increased via a refinement process
	! 2. In the hard case, which has a solution with a special form
	! 3. Error in above search, terminate with failure

	if (found_lambdaC_in_L)
	{
		! Refine an estimate lambdaC in L
		Number lambda_plus, lambda1_neg1, lambda3_2;
		DFOPT_LOG_VERBOSE_IF(debug) << "Phase 2: refining lambdaC = " << lambdaC;
		for (Index it2 = 0; it2 < num_refinement_iters; ++it2)
		{
			status = taylor_new_lambda(x, lambdaC, delta, -1.0, 1, lambda1_neg1, is_tr);
			if (status == 0)
			{
				status = taylor_new_lambda(x, lambdaC, delta, 2.0, 3, lambda3_2, is_tr);
				!DFOPT_LOG_VERBOSE_IF(debug) << "It " << it << ", lambda1(-1) = " << lambda1_neg1 << ", lambda3(2) = " << lambda3_2;
				if (status == 0)
					lambda_plus = std::max(lambda1_neg1, lambda3_2);
				else
					lambda_plus = lambda1_neg1;
			}
			else
			{
				lambda_plus = lambdaC;
			}
			DFOPT_LOG_VERBOSE_IF(debug) << "- Refining iteration " << it << " found new lambdaC <-" << lambda_plus;
			if (std::abs(lambdaC - lambda_plus) < EPS_MACHINE * std::max(1.0, std::abs(lambdaC)))
			{
				! termination from GALAHAD/trs.f90 -- refinement iteration not achieveing much
				lambdaC = lambda_plus;
				break;
			}
			lambdaC = lambda_plus;

			! Recompute factorization
			status = solve_kkt_system(H_scal, g_scal, lambdaC, x, dval);
			if (status == 0)
			{
				! H + lambdaC*I is positive definite
				! Check if we are near the desired norm
				norm_check = (is_tr ? delta : lambdaC / delta);  ! check ||x(lambda)|| ~ this value
				if (std::abs(x.norm2() - norm_check) < boundary_thresh * std::max(1.0, norm_check))
				{
					DFOPT_LOG_VERBOSE_IF(debug) << "Terminating near boundary (success)";
					break;
				}
			}
			else
			{
				! In the refinement phase (lambdaC in L), we should never be able to produce a new lambdaC in N
				! i.e. H + lambdaC*I should always be positive definite
				if (status > 0)
					! H + lambdaC*I is not positive definite
					DFOPT_LOG_VERBOSE_IF(debug) << "ISSUE, lambdaC has left the good region!";
					else
						! Error in solve_kkt_system()
					DFOPT_LOG_VERBOSE_IF(debug) << "Positive definite Cholesky solve failed, stopping";
				break;
			}
		}

		! End of refinement phase, terminate with good estimate...
		if (status == 0)
		{
			DFOPT_LOG_VERBOSE_IF(debug) << "Final solve with lambdaC = " << lambdaC;
			status = solve_kkt_system(H_scal, g_scal, lambdaC, x, dval);
			if (status == 0)
			{
				Number xnorm = x.norm2();
				norm_check = (is_tr ? delta : lambdaC / delta);  ! check ||x(lambda)|| ~ this value
				if (std::abs(xnorm - norm_check) < boundary_thresh * std::max(1.0, norm_check))
				{
					DFOPT_LOG_VERBOSE_IF(debug) << "Phase 2 success";
					result = true;
				}
				else
				{
					! Scale x to have the desired norm
					DFOPT_LOG_VERBOSE_IF(debug) << "Phase 2 solution too far from boundary, ||x|| - expected value = " << xnorm - norm_check << ", scaling to desired norm";
					x.scale(norm_check / xnorm);  ! x <-- (norm_check / xnorm)
					result = true;
				}
			}
			else
			{
				if (status > 0)
				{
					DFOPT_LOG_VERBOSE_IF(debug) << "Easy case error: final lambdaC gave indefinite Hessian";
				}
				else
				{
					DFOPT_LOG_VERBOSE_IF(debug) << "Positive definite Cholesky solve failed, stopping";
				}
				result = false;
			}
		}
		else
		{
			! Error in refinement phase
			DFOPT_LOG_VERBOSE_IF(debug) << "Error in refinement phase";
			result = false;
		}
	}
	else if (hard_case)
	{
		DFOPT_LOG_VERBOSE_IF(debug) << "Hard case";
		! Here, z is a good estimate of u1, a unit eigenvector corresponding to lambda1
		lambdaC = -H_scal.quadform(z);  ! lambdaC = -lambda1
		status = solve_kkt_system(H_scal, g_scal, lambdaC, x, dval);

		if (status > 0)
		{
			! H + lambdaC*I not positive definite, try again with slightly larger lambdaC to avoid rounding errors
			lambdaC += std::max(rounding_perturbation, rounding_perturbation * std::abs(lambdaC));
			status = solve_kkt_system(H_scal, g_scal, lambdaC, x, dval);
		}

		if (status == 0)
		{
			! Final solution is x + alpha*z, with alpha chosen to give correct vector norm
			Number alpha;
			status = hard_case_stepsize(x, z, delta, lambdaC, alpha, is_tr);
			if (status == 0)
			{
				DFOPT_LOG_VERBOSE_IF(debug) << " - Using alpha = " << alpha;
				x.add_multiple(alpha, z);  ! x += alpha*z
			}
			! status != 0, i.e. if no roots to the quadratic, then ||xs|| sufficiently large already, so nothing to do
			result = true;
		}
		else
		{
			if (status > 0)
			{
				DFOPT_LOG_VERBOSE_IF(debug) << "FAILURE: Rounding errors, Rayleigh quotient gave eigenvalue underestimate";
			}
			else
			{
				DFOPT_LOG_VERBOSE_IF(debug) << "Positive definite Cholesky solve failed, stopping";
			}
			result = false;
		}
	}  ! end hard case
	else
	{
		DFOPT_LOG_VERBOSE_IF(debug) << "Error in bisection phase";
		result = false;
	}

	! Catch any inf/NaN issues here
	if (!x.all_finite()) result = false;

	! Compute Cauchy step as safeguard - use if above computation failed, or didn't get sufficient decrease
	cauchy_step(H_scal, g_scal, delta, is_tr);  ! set xcauchy to Cauchy step
	bool use_cauchy_step = false;

	if (result)
	{
		! If global minimizer computation above succeeded, make sure got at least Cauchy decrease
		use_cauchy_step = false;
		Number cauchy_decrease = (is_tr ? model_decrease(H_scal, g_scal, xcauchy) : model_decrease_cubic_reg(H_scal, g_scal, delta, xcauchy));
		Number current_decrease = (is_tr ? model_decrease(H_scal, g_scal, x) : model_decrease_cubic_reg(H_scal, g_scal, delta, x));
		
		if (cauchy_decrease > current_decrease)
		{
			DFOPT_LOG_VERBOSE_IF(debug) << "Global step didn't achieve sufficient decrease, using Cauchy step instead";
			use_cauchy_step = true;
		}
		else if (is_tr && (x.norm2() > delta + boundary_thresh * std::max(1.0, delta)))
		{
			DFOPT_LOG_VERBOSE_IF(debug) << "Global step outside feasible region, using Cauchy step instead";
			use_cauchy_step = true;
		}
		else
		{
			use_cauchy_step = false;
		}
	}
	else
	{
		DFOPT_LOG_VERBOSE_IF(debug) << "Global step calculation failed, using Cauchy step instead";
		use_cauchy_step = true;
	}

	if (use_cauchy_step) x.copy_from(xcauchy);  ! x <-- copy(xcauchy)

	return true;
end subroutine trglob


subroutine initlda(const SymmetricMatrix& H, const Vector& g, const Number delta,
	Number& lambdaL, Number& lambdaU, const bool is_tr)
	const Number expand_interval_thresh = 1e-5; ! expand [lambdaL, lambdaU] interval slightly, in case correct lambda is at endpoint

	! Bounds on min/max eigenvalues of H, satisfying:
	! -lambda_min(H) <= lambda_min_bound <-- NOTE SIGN FLIP
	! lambda_max(H) <= lambda_max_bound
	Number lambda_min_bound, lambda_max_bound;

	! Use Frobenius and infinity norms of H to get first estimates on min/max eigenvalues
	Number H_normF = H.normF();
	Number H_normInf = H.normInf();
	Number H_min_diag = H.get(0, 0);
	for (Index i = 1; i < H.dim(); ++i)  ! start from i=1
		H_min_diag = std::min(H_min_diag, H.get(i, i));

	lambda_min_bound = std::min(H_normF, H_normInf);
	lambda_max_bound = std::min(H_normF, H_normInf);

	! Use Gershgorin discs to get next estimates on min/max eigenvalues
	Number gershgorin_lower_bound, gershgorin_upper_bound; ! gershgorin_lower_bound <= lambda(H) <= gershgorin_upper_bound
	H.eigenvalue_bounds(gershgorin_lower_bound, gershgorin_upper_bound);
	lambda_min_bound = std::min(lambda_min_bound, -gershgorin_lower_bound);  ! note sign flip
	lambda_max_bound = std::min(lambda_max_bound, gershgorin_upper_bound);

	! Norm of g
	Number normg = g.norm2();

	! Set bounds
	lambdaL = std::max(0.0, -H_min_diag);
	lambdaU = 0.0;
	if (is_tr)
	{
		lambdaL = std::max(lambdaL, normg / delta - lambda_max_bound);
		lambdaU = std::max(lambdaU, normg / delta + lambda_min_bound);
	}
	else
	{
		! Bounds given by largest root of: lambda^2 + lambda_{min/max} * lambda - ||g|| * delta = 0
		! Always has at exactly one positive root when ||g||*delta > 0
		! Always have delta>0 from earlier checks, so if ||g||=0 then either H convex -> lambdaC=0 -> x=0
		! or H indefinite -> hard case lambdaC=-lambda1
		! Recall lambda_min_bound = -lambda_min estimate
		if (normg == 0.0)
		{
			lambdaL = std::max(lambdaL, lambda_min_bound); ! lambdaL = -ve best estimate of lambda_min, or zero
			lambdaU = lambdaL + std::max(expand_interval_thresh, expand_interval_thresh * lambdaL);  ! lambdaU slightly above lambdaL
		}
		else
		{
			Index nroots;
			Number l1, l2;
			nroots = quadroots(1.0, -lambda_min_bound, -normg * delta, l1, l2);
			!std::cout << "nroots = " << nroots << ", vals = " << l1 << ", " << l2 << std::endl;
			if (nroots == 2) lambdaU = std::max(lambdaU, std::max(l1, l2));

			nroots = quadroots(1.0, lambda_max_bound, -normg * delta, l1, l2);
			!std::cout << "nroots = " << nroots << ", vals = " << l1 << ", " << l2 << std::endl;
			if (nroots == 2) lambdaL = std::max(lambdaL, std::max(l1, l2));
		}
	}

	! Just in case the true lambdaC is exactly equal to one of these, widen the interval slightly so it becomes an interior point
	if (lambdaL > 0) lambdaL -= std::max(expand_interval_thresh, expand_interval_thresh * lambdaL);
	lambdaL = std::max(lambdaL, 0.0);  ! the above may make lambdaL slightly negative, which we don't want
	if (lambdaU > 0) lambdaU += std::max(expand_interval_thresh, expand_interval_thresh * lambdaU);
end subroutine initlda

subroutine solvekkt(const SymmetricMatrix& H, const Vector& g, 
	const Number lambdaC, Vector& x, Number& dval)
	H_plus_lambda_I.copy_from(H);  ! H_plus_lambda_I <-- copy(H)
	H_plus_lambda_I.add_multiple_identity(lambdaC);  ! H_plus_lambda_I += lambdaC*I
	Index status = H_plus_lambda_I.cholesky_factorize_inplace_safe(x, dval);
	if (status == 0)
	{
		! H + lambdaC*I was positive definite: solve (H + lambdaC*I) * x = -g
		x.copy_from(g);  ! x <-- copy(g)
		x.scale(-1.0);  ! x *= -1
		H_plus_lambda_I.cholesky_solve(x);
	}
	return status;
end subroutine solvekkt

subroutine pi3(const Vector& x_lambda, Number& pi, 
	Number& d1pi, Number& d2pi, Number& d3pi)
	pi = x_lambda.sqnorm2(); ! zero-th deriv is just pi(lambda) = ||x(lambda)||^2

	const Number alpha0 = 1.0;
	const Number alpha1 = 6.0;

	x1.copy_from(x_lambda);  ! x1 <-- copy(x_lambda)
	x1.scale(-1.0);  ! x1 *= -1
	H_plus_lambda_I.cholesky_solve(x1);  ! x1 <-- (H+lambda*I) \ x1 = -1 * (H+lambda*I) \ x_lambda
	x2.copy_from(x1);  ! x2 <-- copy(x1)
	x2.scale(-2.0);  ! x2 *= -2
	H_plus_lambda_I.cholesky_solve(x2);  ! x2 <-- (H+lambda*I) \ x2 = -2 * (H+lambda*I) \ x1

	d1pi = 2.0 * alpha0 * dot(x_lambda, x1); ! pi'(lambda) = 2*alpha0 * dot(x_lambda, x1)
	d2pi = alpha1 * x1.sqnorm2(); ! pi''(lambda) = alpha1 * ||x1||^2
	d3pi = 2.0 * alpha1 * dot(x1, x2); ! pi'''(lambda) = 2*alpha1 * dot(x1, x2)
endsubroutine pi3

subroutine pi3beta(const Vector& x_lambda, const Number beta, 
	Number& pi_beta, Number& d1pi_beta, Number& d2pi_beta, Number& d3pi_beta)
	Number pi, d1pi, d2pi, d3pi;
	Number half_beta = 0.5 * beta;
	pi_three_derivs(x_lambda, pi, d1pi, d2pi, d3pi);
	pi_beta = std::pow(pi, half_beta);
	d1pi_beta = half_beta * std::pow(pi, half_beta - 1) * d1pi;
	d2pi_beta = half_beta * std::pow(pi, half_beta - 1) * d2pi
		+ half_beta * (half_beta - 1) * std::pow(pi, half_beta - 2) * std::pow(d1pi, 2);
	d3pi_beta = std::pow(pi, 2) * d3pi + 3 * (half_beta - 1) * pi * d1pi * d2pi
		+ (half_beta - 1) * (half_beta - 2) * std::pow(d1pi, 3);
	d3pi_beta *= half_beta * std::pow(pi, half_beta - 3);
end subroutine pi3beta

subroutine newlda(const Vector& x_lambda, const Number lambdaC, 
	const Number delta, const Number beta, const Index k, Number& new_lambda, const bool is_tr)
	! Check for valid combinations of (beta,k)
	if (k < 1 || k > 3 || beta == 0.0) return -1;
	if (!is_tr)
		if (!(beta == -1.0 && k == 1) && !(beta == 2.0 && k == 2) && !(beta == 2.0 && k == 3)) return -1;

	Number pi_beta, d1pi_beta, d2pi_beta, d3pi_beta;
	Number sq_delta = delta * delta;
	pi_three_derivs_with_beta(x_lambda, beta, pi_beta, d1pi_beta, d2pi_beta, d3pi_beta);

	Number d1, d2, d3;
	Index nroots;

	if (k == 1)
	{
		! For (TRS), solve: pi_beta + d1pi_beta*d - delta^beta = 0
		!
		! For (ARC), then have beta=-1, and solve
		!    pi_beta + d1pi_beta*d - (lambdaC+d)^(-1) / delta^(-1) = 0
		! or
		!    d1pi_beta * d^2 + (pi_beta + d1pi_beta*lambdaC) * d + (lambdaC*pi_beta - delta) = 0
		if (is_tr)
			nroots = cubic_roots(0.0, 0.0, d1pi_beta, pi_beta - std::pow(delta, beta), d1, d2, d3);
		else
			nroots = cubic_roots(0.0, d1pi_beta, pi_beta + d1pi_beta * lambdaC, lambdaC * pi_beta - delta, d1, d2, d3);
	}
	else if (k == 2)
	{
		! For (TRS), solve: pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 - delta^beta = 0
		!
		! For (ARC), then have beta=2, so solve:
		!   pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 - (lambdaC+d)^2 / delta^2 = 0
		! or
		!   (pi_beta - lambdaC^2 / delta^2) + (d1pi_beta - 2*lambdaC/delta^2) * d + (0.5*d2pi_beta - 1/delta^2) * d^2 = 0
		if (is_tr)
			nroots = cubic_roots(0.0, 0.5 * d2pi_beta, d1pi_beta, pi_beta - std::pow(delta, beta), d1, d2, d3);
		else
			nroots = cubic_roots(0.0, 0.5 * d2pi_beta - 1.0 / (sq_delta), d1pi_beta - 2.0 * lambdaC / sq_delta, pi_beta - lambdaC * lambdaC / sq_delta, d1, d2, d3);
	}
	else
	{
		! For (TRS), solve: pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 + (1/6) * d3pi_beta * d^3 - delta^beta = 0
		!
		! For (ARC), then have beta=2, so solve:
		!     pi_beta + d1pi_beta*d + 0.5*d2pi_beta * d^2 + (1/6) * d3pi_beta * d^3 - (lambdaC+d)^2 / delta^2 = 0
		! or
		!     (pi_beta - lambdaC^2 / delta^2) + (d1pi_beta - 2*lambdaC/delta^2) * d + (0.5*d2pi_beta - 1/delta^2) * d^2 + (1/6) * d3pi_beta * d^3 = 0
		if (is_tr)
			nroots = cubic_roots(d3pi_beta / 6.0, 0.5 * d2pi_beta, d1pi_beta, pi_beta - std::pow(delta, beta), d1, d2, d3);
		else
			nroots = cubic_roots(d3pi_beta / 6.0, 0.5 * d2pi_beta - 1.0 / (sq_delta), d1pi_beta - 2.0 * lambdaC / sq_delta, pi_beta - lambdaC * lambdaC / sq_delta, d1, d2, d3);
	}

	if (nroots == 0) {
		new_lambda = lambdaC;
		return -1;
	}
	else if (nroots == 1)
		new_lambda = lambdaC + d1;
	else if (nroots == 2)
		new_lambda = lambdaC + d2;
	else
		new_lambda = lambdaC + d3;

	return 0;
end subroutine newlda

subroutine hardstep(const Vector& xs, const Vector& ztmp, 
	const Number delta, const Number lambdaC, Number& alpha, const bool is_tr)
	Number alpha1, alpha2;
	Number rhs = (is_tr ? delta : lambdaC / delta);
	Index nroots = quadroots(ztmp.sqnorm2(), 2.0 * dot(xs, ztmp), xs.sqnorm2() - rhs * rhs, alpha1, alpha2);
	if (nroots > 0)
	{
		alpha = alpha1;  ! any root is fine, since sign of z is arbitrary
		return 0;
	}
	else
		return -1;
end subroutine hardstep

subroutine cauchy(hess, g, delta, is_tr, s)
    ! Common modules
    use, non_intrinsic :: consts_mod, only : RP, IK, HALF, DEBUGGING
    use, non_intrinsic :: debug_mod, only : assert
    use, non_intrinsic :: linalg_mod, only : matprod, inprod, issymmetric

    implicit none

    ! Inputs
    real(RP), intent(in):: hess(:, :)  ! HESS(N,N)
    real(RP), intent(in):: g(:)  ! G(N,N)
    real(RP), intent(in):: delta
    logical, intent(in) :: is_tr

    ! Outputs
    real(RP), intent(out) :: s

    ! Local variables
    character(len=*), parameter :: srname = 'CAUCHY'
    integer(IK) :: n
    real(RP) :: alpha
    real(RP) :: gHg
    real(RP) :: normg
    real(RP) :: tmp
    real(RP) :: hg(size(g))

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
    int(IP), intent(out) :: nroots

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
        if (p2 == ZERO)
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
		if (dp .neq. ZERO) then
            x1 -= p / dp;
        end if
    end if

	if (nroots >= 2) then
		! refine x2
		p = p0 * x2 * x2 + p1 * x2 + p2;
		dp = 2 * p0 * x2 + p1;
		if (dp .neq. ZERO) then
            x2 -= p / dp;
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
    use, non_intrinsic :: consts_mod, only : RP, IK, HALF, ZERO

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
    int(IP), intent(out) :: nroots

    ! Locals
    real(RP) :: q1, q2, q3
    real(RP) :: s, t, b, c, d
    real(RP) :: p, dp

	if (p0 == ZERO) then
		call quadroots(p1, p2, p3, x1, x2, nroots)
	else if (p3 == ZERO)
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
			if (d .neq. ZERO) then
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
		if (dp .neq. ZERO) then
            x1 -= p / dp
        end if
    end if

	if (nroots >= 2) then
		! refine x2
		p = p0 * x2 * x2 * x2 + p1 * x2 * x2 + p2 * x2 + p3
		dp = 3 * p0 * x2 * x2 + 2 * p1 * x2 + p2
		if (dp .neq. ZERO) then
            x2 -= p / dp
        end if
    end if

	if (nroots >= 3) then
		! refine x3
		p = p0 * x3 * x3 * x3 + p1 * x3 * x3 + p2 * x3 + p3
		dp = 3 * p0 * x3 * x3 + 2 * p1 * x3 + p2
		if (dp .neq. ZERO) then
            x3 -= p / dp
        end if
    end if

end subroutine cubicroots

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