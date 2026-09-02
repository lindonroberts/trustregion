! Global module to expose main routines

module trustregion_mod

    use, non_intrinsic :: consts_mod, only : RP, IK

    implicit none

    interface
        module subroutine trsapp(delta, g_in, hess_in, tol, crvmin, s, info)
            real(RP), intent(in) :: delta
            real(RP), intent(in) :: g_in(:)   ! G_IN(N)
            real(RP), intent(in) :: hess_in(:, :)    ! HESS_IN(N, N)
            real(RP), intent(in) :: tol
            real(RP), intent(out) :: crvmin
            real(RP), intent(out) :: s(:)   ! S(N)
            integer(IK), intent(out), optional :: info
        end subroutine trsapp

        module subroutine trsbox(delta, g_in, hess_in, sl, su, tol, xopt, crvmin, d)
            real(RP), intent(in) :: delta
            real(RP), intent(in) :: g_in(:)  ! G_IN(N)
            real(RP), intent(in) :: hess_in(:, :)  ! HESS_IN(N, N)
            real(RP), intent(in) :: sl(:)  ! SL(N)
            real(RP), intent(in) :: su(:)  ! SU(N)
            real(RP), intent(in) :: tol
            real(RP), intent(in) :: xopt(:)  ! XOPT(N)
            real(RP), intent(out) :: crvmin
            real(RP), intent(out) :: d(:)  ! D(N)
        end subroutine trsbox

        module subroutine trslin(amat_in, bvec_in, xopt, delta, g_in, hess_in, tol, s)
            real(RP), intent(in) :: amat_in(:, :)  ! AMAT_IN(N, M)
            real(RP), intent(in) :: bvec_in(:)  ! BVEC_IN(M)
            real(RP), intent(in) :: xopt(:)  ! XOPT(N)
            real(RP), intent(in) :: delta
            real(RP), intent(in) :: g_in(:)  ! G_IN(N)
            real(RP), intent(in) :: hess_in(:, :)  ! HESS_IN(N, N)
            real(RP), intent(in) :: tol
            real(RP), intent(out) :: s(:)  ! S(N)
        end subroutine trslin
    end interface
end module trustregion_mod