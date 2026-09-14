module trustregion_c_api

    use, intrinsic :: iso_c_binding
    use, non_intrinsic :: trustregion_mod, only : trsunc, arcunc, trsapp, trsbox, trslin

    implicit none

contains

    subroutine c_trsunc(n, delta, g_in, hess_in, lambda, s) & 
        bind(C, name="f_trsunc")
        integer(c_int), value :: n
        real(c_double), intent(in) :: delta
        real(c_double), intent(in) :: g_in(n)
        real(c_double), intent(in) :: hess_in(n * n)
        real(c_double), intent(out) :: lambda
        real(c_double), intent(out) :: s(n)

        ! Locals
        integer(c_int) :: i, j
        real(c_double) :: hess(n, n)

        ! Convert c-style array hess_in to Fortran matrix
        do i = 1, n 
            do j = 1, n 
                hess(i, j) = hess_in((i-1) * n + j)
            end do
        end do

        call trsunc(delta, g_in, hess, lambda, s)
    end subroutine c_trsunc

    subroutine c_arcunc(n, delta, g_in, hess_in, lambda, s) & 
        bind(C, name="f_arcunc")
        integer(c_int), value :: n
        real(c_double), intent(in) :: delta
        real(c_double), intent(in) :: g_in(n)
        real(c_double), intent(in) :: hess_in(n * n)
        real(c_double), intent(out) :: lambda
        real(c_double), intent(out) :: s(n)

        ! Locals
        integer(c_int) :: i, j
        real(c_double) :: hess(n, n)

        ! Convert c-style array hess_in to Fortran matrix
        do i = 1, n 
            do j = 1, n 
                hess(i, j) = hess_in((i-1) * n + j)
            end do
        end do

        call arcunc(delta, g_in, hess, lambda, s)
    end subroutine c_arcunc

    subroutine c_trsapp(n, delta, g_in, hess_in, tol, crvmin, s, info) & 
        bind(C, name="f_trsapp")
        integer(c_int), value :: n
        real(c_double), intent(in) :: delta
        real(c_double), intent(in) :: g_in(n)
        real(c_double), intent(in) :: hess_in(n * n)
        real(c_double), intent(in) :: tol
        real(c_double), intent(out) :: crvmin
        real(c_double), intent(out) :: s(n)
        integer(c_int), intent(out) :: info

        ! Locals
        integer(c_int) :: i, j
        real(c_double) :: hess(n, n)

        ! Convert c-style array hess_in to Fortran matrix
        do i = 1, n 
            do j = 1, n 
                hess(i, j) = hess_in((i-1) * n + j)
            end do
        end do

        !print *, "calling trsapp..."
        !print *, "delta = ", delta
        !print *, "g_in = ", g_in
        !print *, "hess_in = ", hess_in
        !print *, "tol = ", tol

        !print *, "hess ="
        !do i = 1, n 
        !    print *, hess(i, :)
        !end do

        call trsapp(delta, g_in, hess, tol, crvmin, s, info)

        !print *, "after trsapp..."
        !print *, "crvmin = ", crvmin
        !print *, "s = ", s 
        !print *, "info = ", info
    end subroutine c_trsapp

    subroutine c_trsbox(n, delta, g_in, hess_in, sl, su, tol, xopt, crvmin, s) & 
        bind(C, name="f_trsbox")
        integer(c_int), value :: n
        real(c_double), intent(in) :: delta
        real(c_double), intent(in) :: g_in(n)
        real(c_double), intent(in) :: hess_in(n * n)
        real(c_double), intent(in) :: sl(n)
        real(c_double), intent(in) :: su(n)
        real(c_double), intent(in) :: tol
        real(c_double), intent(in) :: xopt(n)
        real(c_double), intent(out) :: crvmin
        real(c_double), intent(out) :: s(n)

        ! Locals
        integer(c_int) :: i, j
        real(c_double) :: hess(n, n)

        ! Convert c-style array hess_in to Fortran matrix
        do i = 1, n 
            do j = 1, n 
                hess(i, j) = hess_in((i-1) * n + j)
            end do
        end do

        call trsbox(delta, g_in, hess, sl, su, tol, xopt, crvmin, s)
    end subroutine c_trsbox

    subroutine c_trslin(n, m, amat_in, bvec_in, xopt, delta, g_in, hess_in, tol, s) & 
        bind(C, name="f_trslin")
        integer(c_int), value :: n
        integer(c_int), value :: m
        real(c_double), intent(in) :: amat_in(n * m)
        real(c_double), intent(in) :: bvec_in(m)
        real(c_double), intent(in) :: xopt(n)
        real(c_double), intent(in) :: delta
        real(c_double), intent(in) :: g_in(n)
        real(c_double), intent(in) :: hess_in(n * n)
        real(c_double), intent(in) :: tol
        real(c_double), intent(out) :: s(n)

        ! Locals
        integer(c_int) :: i, j
        real(c_double) :: hess(n, n)
        real(c_double) :: amat(n, m)

        ! Convert c-style array hess_in to Fortran matrix
        do i = 1, n 
            do j = 1, n 
                hess(i, j) = hess_in((i-1) * n + j)
            end do
        end do

        do i = 1, n 
            do j = 1, m
                amat(i, j) = amat_in((i-1) * n + j)
            end do
        end do

        call trslin(amat, bvec_in, xopt, delta, g_in, hess, tol, s)
    end subroutine c_trslin


end module trustregion_c_api