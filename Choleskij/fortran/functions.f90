subroutine Choleskij_Decomposition(A, n, comp_mode, b, L, time, error)
        integer(4) :: n, comp_mode, b, i, info, m
        real(8), dimension(n,n) :: A, L
        integer(4), allocatable, dimension(:) :: ipiv
        real(8), allocatable, dimension(:) :: work
        real(8), allocatable, dimension(:) :: L_res
        real(8), allocatable, dimension(:,:) :: Lt, E
        real(8) :: time, error, start, finish
        double precision zero, one
        parameter (zero = 0.0D0, one = 1.0D0)

        L = 0

        call cpu_time(start)

        !Fortran default
        if (comp_mode == 0) then
                do i = 1, n

                end do

        !BLAS
        else if (comp_mode == 1) then
                allocate(L_res(n))
                L_res = 0
                do i = 1, n

                end do
        deallocate(L_res)

        !LAPACK
        else if (comp_mode == 2) then
                L = A
                allocate(ipiv(n))
                !call dpbtrf(n, n, LU, n, ipiv, info)
        end if

        call cpu_time(finish)
        time = finish - start

        if (comp_mode == 2) then !||A*A^(-1) - I||
                allocate(work(n*n))
                allocate(E(n,n))
                do i = 1, n
                        E(i,i) = 1
                end do
                call dgetri(n, L, n, ipiv, work, n, info)
                error = norm2(matmul(A, L) - E)
                deallocate(ipiv)
                deallocate(work)
                deallocate(E)
        else !||LU - A||
                        error = norm2(matmul(L, Lt) - A)
    end if
end subroutine
