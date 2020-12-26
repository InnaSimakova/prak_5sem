        program LLtdecomposition
        use header
        integer(4) :: n, comp_mode, b
        character(8) :: manual_mode
        real(8) :: time, error
        real(8), allocatable, dimension(:,:) :: A, L

        print '(a)', "Enter square matrix (n x n) dimension:"
        print '(a,$)', "n = "
        read *, n
        allocate(A(n,n))
        allocate(L(n,n))
        L = 0
        b = 0

        do
                print '(a)', "Manual input of matrix? (y/n)"
                read *, manual_mode

        if (manual_mode == 'y') then
            print "(a, i5, a)", "Enter square ", n !, "-dimensional
            !matrix:"
            read *, A
            A = transpose(A)
            exit
        else if (manual_mode == 'n') then
                call random_number(A)
                A = floor(20 * A)
                A = A + 1
                exit
                else
                        print '(a)', "Answer not recognized, try again."
                end if
        end do


        do
                print '(a)', "Choose computation: (type 0,1 or 2)"
                print '(a)', "0 - default Fortran functions."
                print '(a)', "1 - BLAS library functions."
                print '(a)', "2 - LAPACK library functions."
                read *, comp_mode

                if ((comp_mode == 0) .or. (comp_mode == 1) .or. (comp_mode == 2) .or. (comp_mode == 3)) then
                        exit
                else
                        print '(a)', "Error! Try again."
                end if
        end do

        print "(6f10.3)", A

        call Choleskij_Decomposition(A, n, comp_mode, b, L, time, error)

        print "(6f10.3)", transpose(L)
        print '("Time = ", f10.2," ms")', 1000 * time
        print '("||LU - A|| = ", E10.2)', error

        deallocate(A)
        deallocate(L)
        end program
