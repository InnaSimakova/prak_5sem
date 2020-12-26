        module header
                interface
                subroutine Choleskij_Decomposition(A, n, comp_mode, b, L, time, error)
                        integer(4) :: n, comp_mode, b
                        real(8), dimension(n,n) :: A, L
                        real(8) :: time, error
                end subroutine
                end interface
        end module
